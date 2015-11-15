#!/usr/bin/perl -w

use English;
use warnings;
use strict;
use Carp qw(carp cluck croak confess);
use Getopt::Long qw( :config posix_default bundling no_ignore_case );
use POSIX qw(ceil floor);
use List::Util qw[min max];
use Cwd 'abs_path';
use Cwd;

use cworld::dekker;

my $tool=(split(/\//,abs_path($0)))[-1];

sub check_options {
    my $opts = shift;

    my ($inputMatrix,$mergeSelf,$mergeMode,$verbose,$output,$excludeZero,$fillGaps);
    
    my $ret={};
    
    if( exists($opts->{ inputMatrix }) ) {
        $inputMatrix = $opts->{ inputMatrix };
    } else {
        print STDERR "\nERROR: Option inputMatrix|i is required.\n";
        help();
    }
   
   if( exists($opts->{ mergeSelf }) ) {
        $mergeSelf = $opts->{ mergeSelf };
    } else {
        $mergeSelf=undef
    }
    
    if( exists($opts->{ mergeMode }) ) {
        $mergeMode = $opts->{ mergeMode };
    } else {
        $mergeMode="mean";
    }
    
    if( exists($opts->{ verbose }) ) {
        $verbose = $opts->{ verbose };
    } else {
        $verbose = 0;
    }
    
    if( exists($opts->{ output }) ) {
        $output = $opts->{ output };
    } else {
        $output = "";
    }
    
    if( exists($opts->{ excludeZero }) ) {
        $excludeZero = 1;
    } else {
        $excludeZero = 0;
    }
    
    if( exists($opts->{ fillGaps }) ) {
        $fillGaps = 1;
    } else {
        $fillGaps = 0;
    }
    
    $ret->{ inputMatrix }=$inputMatrix;
    $ret->{ mergeSelf }=$mergeSelf;
    $ret->{ mergeMode }=$mergeMode;
    $ret->{ verbose }=$verbose;
    $ret->{ output }=$output;
    $ret->{ excludeZero }=$output;
    $ret->{ fillGaps }=$fillGaps;
    
    return($ret,$inputMatrix,$mergeSelf,$mergeMode,$verbose,$output,$excludeZero,$fillGaps);
}


sub intro() {
    print STDERR "\n";
    
    print STDERR "Tool:\t\t".$tool."\n";
    print STDERR "Version:\t".$cworld::dekker::VERSION."\n";
    print STDERR "Summary:\ttransform rectangular matrix into square (symmetrical) matrix\n";
    
    print STDERR "\n";
}

sub help() {
    intro();
    
    print STDERR "Usage: perl matrix2symmetrical.pl [OPTIONS] -i <inputMatrix>\n";
    
    print STDERR "\n";
    
    print STDERR "Required:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-i", "[]", "input matrix file");
    
    print STDERR "\n";
    
    print STDERR "Options:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-v", "[]", "FLAG, verbose mode");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-o", "[]", "prefix for output file(s)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ms", "[]", "FLAG, merge self fragment interactions");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--mm", "[]", "merge mode, method to merge signal (mean,median,min,max,etc)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ez", "[]", "FLAG, ignore 0s in all calculations");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--fg", "[]", "FLAG, add row/cols for missing primers/fragments (gaps)");
    
    print STDERR "\n";
    
    print STDERR "Notes:";
    print STDERR "
    This script turns any matrix into a symmetrical matrix.  If the matrix is already symmetrical, nothing is changed.
    Matrix can be TXT or gzipped TXT.
    See website for matrix format details.\n";
    
    print STDERR "\n";
    
    print STDERR "Contact:
    Bryan R. Lajoie
    Dekker Lab 2015
    https://github.com/blajoie/cworld-dekker
    http://my5C.umassmed.edu";
    
    print STDERR "\n";
    print STDERR "\n";
    
    exit;
}

sub symmetricizeData($$) {
    my $matrixObject=shift;
    my $matrix=shift;
    
    my $inc2header=$matrixObject->{ inc2header };
    my $header2inc=$matrixObject->{ header2inc };
    my $numYHeaders=$matrixObject->{ numYHeaders };
    my $numXHeaders=$matrixObject->{ numXHeaders };    
    my $missingValue=$matrixObject->{ missingValue };
    my $verbose=$matrixObject->{ verbose };
    
    print STDERR "\tsymmetrical\n" if($verbose);
    
    for(my $y=0;$y<$numYHeaders;$y++) {
        my $yHeader=$inc2header->{ y }->{$y};
        for(my $x=0;$x<$numXHeaders;$x++) {
            my $xHeader=$inc2header->{ x }->{$x};    
            
            next if($x < $y);

            my $cScore=$matrixObject->{ missingValue };
            die("error, matrix appears symmetrical ($matrix->{$y}->{$x} vs $matrix->{$x}->{$x})!") if((defined($matrix->{$y}->{$x})) and (defined($matrix->{$x}->{$y})));
            $cScore=$matrix->{$y}->{$x} if(defined($matrix->{$y}->{$x}));
            $cScore=$matrix->{$x}->{$y} if(defined($matrix->{$x}->{$y}));
            
            $matrix->{$x}->{$y}=$cScore;
            $matrix->{$y}->{$x}=$cScore;
            
        }
    }
    
    return($matrix);
}

sub mergeSelf($$$$$$) {
    my $matrixObject=shift;
    my $matrix=shift;
    my $fragmentMapping=shift;
    my $mergeSelf=shift;
    my $mergeMode=shift;
    my $excludeZero=shift;
    
    my $inc2header=$matrixObject->{ inc2header };
    my $header2inc=$matrixObject->{ header2inc };
    my $numYHeaders=$matrixObject->{ numYHeaders };
    my $numXHeaders=$matrixObject->{ numXHeaders };
    my $numTotalHeaders=$matrixObject->{ numTotalHeaders };
    my $missingValue=$matrixObject->{ missingValue };
    my $verbose=$matrixObject->{ verbose };
    
    print STDERR "\tmerge self\n" if($verbose);
    
    # build list of unique fragments, map FOR/REV/LREV/LFOR to frag
    my $headerMap={};
    my $tmp_inc2header={};    
    my $tmp_header2inc={};
    my $numHeaders=0;
    foreach my $fragmentRegion ( sort keys %{$fragmentMapping} ) {
        foreach my $fragmentNumber ( sort {$a<=>$b} keys %{$fragmentMapping->{$fragmentRegion}} ) {
        
            my $chromosome="NA";
            my $start="NA";
            my $end="NA";
            my $assembly="NA";
            my $region="NA";
            my $yPrimerType="NA";
            my $xPrimerType="NA";
                        
            if(defined($fragmentMapping->{$fragmentRegion}->{$fragmentNumber}->{ x })) {
                my $header=$fragmentMapping->{$fragmentRegion}->{$fragmentNumber}->{ x };
                my $headerObject=getHeaderObject($header);
                $headerMap->{ x }->{$numHeaders}=$header;
                
                croak "fragment mis-alignment" if(($chromosome ne "NA") and ($chromosome ne $headerObject->{ chromosome }));
                $chromosome = $headerObject->{ chromosome };
                croak "fragment mis-alignment" if(($start ne "NA") and ($start != $headerObject->{ start }));
                $start = $headerObject->{ start };
                croak "fragment mis-alignment" if(($end ne "NA") and ($end != $headerObject->{ end }));
                $end = $headerObject->{ end };
                $assembly = $headerObject->{ assembly };
                $region = $headerObject->{ region };
                $xPrimerType=$headerObject->{ primerType };
            }
            
            if(defined($fragmentMapping->{$fragmentRegion}->{$fragmentNumber}->{ y })) {
                my $header=$fragmentMapping->{$fragmentRegion}->{$fragmentNumber}->{ y };
                my $headerObject=getHeaderObject($header);
                $headerMap->{ y }->{$numHeaders}=$header;
                
                croak "fragment mis-alignment" if(($chromosome ne "NA") and ($chromosome ne $headerObject->{ chromosome }));
                $chromosome = $headerObject->{ chromosome };
                croak "fragment mis-alignment" if(($start ne "NA") and ($start != $headerObject->{ start }));
                $start = $headerObject->{ start };
                croak "fragment mis-alignment" if(($end ne "NA") and ($end != $headerObject->{ end }));
                $end = $headerObject->{ end };
                $assembly = $headerObject->{ assembly };
                $region = $headerObject->{ region };
                $yPrimerType=$headerObject->{ primerType };
            }
    
            my $primerType="NA";
            $primerType=$yPrimerType if($xPrimerType eq "NA");
            $primerType=$xPrimerType if($yPrimerType eq "NA");
            $primerType=$yPrimerType."-".$xPrimerType if(($yPrimerType ne "NA") and ($xPrimerType ne "NA"));
            $primerType=$yPrimerType if($yPrimerType eq $xPrimerType);
            confess "error! [$region | $fragmentNumber] ($yPrimerType | $xPrimerType)\n" if($primerType eq "NA");
            
            my $header="5C_".$region."_".$primerType."_".$fragmentNumber."|".$assembly."|".$chromosome.":".$start."-".$end;
            $tmp_inc2header->{ x }->{$numHeaders}=$header;
            $tmp_header2inc->{ x }->{$header}=$numHeaders;
            $tmp_inc2header->{ y }->{$numHeaders}=$header;
            $tmp_header2inc->{ y }->{$header}=$numHeaders;
            $tmp_inc2header->{ xy }->{$numHeaders}=$header;
            $tmp_header2inc->{ xy }->{$header}=$numHeaders;
            $numHeaders++;
        }
    }
    
    print STDERR "\tre-building matrix ... \n" if($verbose);
    
    # build new matrix
    my $tmpMatrix={};
    for(my $y=0;$y<$numHeaders;$y++) {
        my $yHeader=$tmp_inc2header->{ y }->{$y};
                
        # get original *FOR for fragment
        my $yHeader_forward="NA";
        $yHeader_forward=$headerMap->{ y }->{$y} if(defined($headerMap->{ y }->{$y}));
        my $yHeader_forward_index=-1;
        $yHeader_forward_index=$header2inc->{ y }->{$yHeader_forward} if(defined($header2inc->{ y }->{$yHeader_forward}));
        
        # get original *REV for fragment
        my $yHeader_reverse="NA";
        $yHeader_reverse=$headerMap->{ x }->{$y} if(defined($headerMap->{ x }->{$y}));
        my $yHeader_reverse_index=-1;
        $yHeader_reverse_index=$header2inc->{ x }->{$yHeader_reverse} if(defined($header2inc->{ x }->{$yHeader_reverse}));
        
        for(my $x=0;$x<$numHeaders;$x++) {
            
            next if($x < $y);
            
            my $xHeader=$tmp_inc2header->{ x }->{$x};    
            
            # get original *FOR for fragment        
            my $xHeader_forward="NA";
            $xHeader_forward=$headerMap->{ y }->{$x} if(defined($headerMap->{ y }->{$x}));
            my $xHeader_forward_index=-1;
            $xHeader_forward_index=$header2inc->{ y }->{$xHeader_forward} if(defined($header2inc->{ y }->{$xHeader_forward}));
        
            # get original *REV for fragment        
            my $xHeader_reverse="NA";
            $xHeader_reverse=$headerMap->{ x }->{$x} if(defined($headerMap->{ x }->{$x}));
            my $xHeader_reverse_index=-1;
            $xHeader_reverse_index=$header2inc->{ x }->{$xHeader_reverse} if(defined($header2inc->{ x }->{$xHeader_reverse}));
            
            
            # extract out original interaction scores
            my $cScore_A="NA";
            if(($yHeader_forward_index != -1) and ($xHeader_reverse_index != -1)) {
                $cScore_A=$matrixObject->{ missingValue };
                $cScore_A=$matrix->{$yHeader_forward_index}->{$xHeader_reverse_index} if(defined($matrix->{$yHeader_forward_index}->{$xHeader_reverse_index}));
            }
            $cScore_A = "NA" if(($cScore_A ne "NA") and (($cScore_A == 0) and ($excludeZero)));
            
            my $cScore_B="NA";
            if(($xHeader_forward_index != -1) and ($yHeader_reverse_index != -1)) {
                $cScore_B=$matrixObject->{ missingValue };
                $cScore_B=$matrix->{$xHeader_forward_index}->{$yHeader_reverse_index} if(defined($matrix->{$xHeader_forward_index}->{$yHeader_reverse_index}));
            }
            $cScore_B = "NA" if(($cScore_B ne "NA") and (($cScore_B == 0) and ($excludeZero)));
            
            my $cScore="NA";
            
            # combine scores
            my @tmpList=();
            push(@tmpList,$cScore_A) if($cScore_A ne "NA");
            push(@tmpList,$cScore_B) if($cScore_B ne "NA");
            
            if(@tmpList > 0) {
                my $tmpListStats=listStats(\@tmpList);
                croak "invalid merge mode - [$mergeMode] not available" if(!defined($tmpListStats->{ $mergeMode }));
                $cScore=$tmpListStats->{ $mergeMode };
            }
            
            # make symmetrical
            $tmpMatrix->{$y}->{$x}=$cScore;
            $tmpMatrix->{$x}->{$y}=$cScore;
            
        }
        my $pcComplete = round((($y/($numHeaders-1))*100),2);
        print STDERR "\e[A" if(($verbose) and ($y != 0));
        printf STDERR "\t%.2f%% complete (".$y."/".($numHeaders-1).")...\n", $pcComplete if($verbose);
    }
    
    $matrixObject->{ inc2header }=$tmp_inc2header;
    $matrixObject->{ header2inc }=$tmp_header2inc;
    $matrixObject->{ numYHeaders }=$numHeaders;
    $matrixObject->{ numXHeaders }=$numHeaders;
    $matrixObject->{ numTotalHeaders }=$numHeaders;
    $matrixObject->{ missingValue }="NA";  
    
    $matrixObject=updateMatrixObject($matrixObject);

    return($tmpMatrix,$matrixObject);
}

sub addMissingPrimers($$) {
    my $matrixObject=shift;
    my $fragmentMapping=shift;
    
    my $inc2header=$matrixObject->{ inc2header };
    my $header2inc=$matrixObject->{ header2inc };
    my $numYHeaders=$matrixObject->{ numYHeaders };
    my $numXHeaders=$matrixObject->{ numXHeaders };
    my $numTotalHeaders=$matrixObject->{ numTotalHeaders };
    my $missingValue=$matrixObject->{ missingValue };
    my $symmetrical=$matrixObject->{ symmetrical };
    
    my $tmp_inc2header={};    
    my $tmp_header2inc={};
    my $numHeaders=0;
    my $gapIndex=0;
    
    for(my $h=0;$h<$numTotalHeaders;$h++) {
        my $header=$inc2header->{ xy }->{$h};
        my $headerObject=getHeaderObject($header);
        my $assembly=$headerObject->{ assembly };
        my $region=$headerObject->{ region };
        my $chromosome=$headerObject->{ chromosome };
        my $start=$headerObject->{ start };
        my $end=$headerObject->{ end };
        my $fragmentNumber=$headerObject->{ fragmentNumber };
            
        $tmp_inc2header->{ x }->{$numHeaders}=$header;
        $tmp_header2inc->{ x }->{$header}=$numHeaders;
        $tmp_inc2header->{ y }->{$numHeaders}=$header;
        $tmp_header2inc->{ y }->{$header}=$numHeaders;
        $tmp_inc2header->{ xy }->{$numHeaders}=$header;
        $tmp_header2inc->{ xy }->{$header}=$numHeaders;
        $numHeaders++;
                
        print "h=$header\n";
        
        next if($h >= ($numTotalHeaders-1));
       
        my $next_header=$inc2header->{ xy }->{($h+1)};
        my $next_headerObject=getHeaderObject($next_header);
        my $next_region=$next_headerObject->{ region };
        my $next_chromosome=$next_headerObject->{ chromosome };
        my $next_start=$next_headerObject->{ start };
        my $next_end=$next_headerObject->{ end };
        my $next_fragmentNumber=$next_headerObject->{ fragmentNumber };
        
        # skip same frag, LFOR/REV | FOR/REV
        next if( ($region eq $next_region) and ($chromosome eq $next_chromosome) and ($start == $next_start) and ($end == $next_end) );

        if( ($region eq $next_region) and ($chromosome eq $next_chromosome) and ($end != $next_start) ) {
            my $gapStart=$end;
            my $gapEnd=$next_start;
            
            my $gapFragmentNumberStart=$fragmentNumber+1;
            my $gapFragmentNumberEnd=$next_fragmentNumber-1;
            my $gapFragmentNumber=(($gapFragmentNumberStart+$gapFragmentNumberEnd)/2);
            $gapFragmentNumber=$gapFragmentNumberStart if($gapFragmentNumberStart == $gapFragmentNumberEnd);
            
            
            my $header="5C_".$region."_GAP_".$gapFragmentNumber."|".$assembly."|".$chromosome.":".$gapStart."-".$gapEnd;
            $gapIndex++;
            
            print "g=$header\n";
            
            $tmp_inc2header->{ x }->{$numHeaders}=$header;
            $tmp_header2inc->{ x }->{$header}=$numHeaders;
            $tmp_inc2header->{ y }->{$numHeaders}=$header;
            $tmp_header2inc->{ y }->{$header}=$numHeaders;
            $tmp_inc2header->{ xy }->{$numHeaders}=$header;
            $tmp_header2inc->{ xy }->{$header}=$numHeaders;
            
            $fragmentMapping->{$region}->{$gapFragmentNumber}->{ x }=$header;
            $fragmentMapping->{$region}->{$gapFragmentNumber}->{ y }=$header;
            
            $numHeaders++;
        }
                
    }
    
    $matrixObject->{ inc2header }=$tmp_inc2header;
    $matrixObject->{ header2inc }=$tmp_header2inc;
    $matrixObject=updateMatrixObject($matrixObject);
    
    return($matrixObject,$fragmentMapping);
}
   
my %options;
my $results = GetOptions( \%options,'inputMatrix|i=s','verbose|v','mergeSelf|ms','mergeMode|mm=s','output|o=s','excludeZero|ez','fillGaps|fg') or croak help();
my ($ret,$inputMatrix,$mergeSelf,$mergeMode,$verbose,$output,$excludeZero,$fillGaps)=check_options( \%options );

intro() if($verbose);

#get the absolute path of the script being used
my $cwd = getcwd();
my $fullScriptPath=abs_path($0);
my @fullScriptPathArr=split(/\//,$fullScriptPath);
@fullScriptPathArr=@fullScriptPathArr[0..@fullScriptPathArr-3];
my $scriptPath=join("/",@fullScriptPathArr);
my $commentLine=getScriptOpts($ret,$tool);

croak "inputMatrix [$inputMatrix] does not exist" if(!(-e $inputMatrix));

# get matrix information
my $matrixObject=getMatrixObject($inputMatrix,$output,$verbose);
my $inc2header=$matrixObject->{ inc2header };
my $header2inc=$matrixObject->{ header2inc };
my $numYHeaders=$matrixObject->{ numYHeaders };
my $numXHeaders=$matrixObject->{ numXHeaders };
my $numTotalHeaders=$matrixObject->{ numTotalHeaders };
my $missingValue=$matrixObject->{ missingValue };
my $symmetrical=$matrixObject->{ symmetrical };
my $equalHeaderFlag=$matrixObject->{ equalHeaderFlag };
my $inputMatrixName=$matrixObject->{ inputMatrixName };
$output=$matrixObject->{ output };

my $fragmentMapping={};
for(my $y=0;$y<$numYHeaders;$y++) {
    my $yHeader=$inc2header->{ y }->{$y};
    my $headerObject=getHeaderObject($yHeader);
    my $region=$headerObject->{ region };
    my $fragmentNumber=$headerObject->{ fragmentNumber };
    $fragmentMapping->{$region}->{$fragmentNumber}->{ y }=$yHeader;
}
for(my $x=0;$x<$numXHeaders;$x++) {
    my $xHeader=$inc2header->{ x }->{$x};  
    my $headerObject=getHeaderObject($xHeader);
    my $region=$headerObject->{ region };
    my $fragmentNumber=$headerObject->{ fragmentNumber };
    $fragmentMapping->{$region}->{$fragmentNumber}->{ x }=$xHeader;
}

# make symmetrical
$inc2header->{ x } = $inc2header->{ xy };
$inc2header->{ y } = $inc2header->{ xy };
$header2inc->{ x } = $header2inc->{ xy };
$header2inc->{ y } = $header2inc->{ xy };
$matrixObject->{ missingValue }="NA";

# update matrix object
$matrixObject->{ inc2header }=$inc2header;
$matrixObject->{ header2inc }=$header2inc;
$matrixObject=updateMatrixObject($matrixObject);

# fill primer/gragment gaps if selected
($matrixObject,$fragmentMapping)=addMissingPrimers($matrixObject,$fragmentMapping) if($fillGaps);
$inc2header=$matrixObject->{ inc2header };
$header2inc=$matrixObject->{ header2inc };
$numYHeaders=$matrixObject->{ numYHeaders };
$numXHeaders=$matrixObject->{ numXHeaders };
$numTotalHeaders=$matrixObject->{ numTotalHeaders };

#read Matrix
my $matrix={};
($matrix)=getData($inputMatrix,$matrixObject,$verbose);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

print STDERR "symmetriciz'n data...\n" if($verbose);
($matrix)=symmetricizeData($matrixObject,$matrix) if(!$mergeSelf);
($matrix,$matrixObject)=mergeSelf($matrixObject,$matrix,$fragmentMapping,$mergeSelf,$mergeMode,$excludeZero) if($mergeSelf);

print STDERR "\n" if($verbose);

my $suffix="symmetrical";
$suffix.=".selfMerged" if($mergeSelf);
$suffix.=".gapped" if($fillGaps);

my $symmetricalMatrixFile = $output.".".$suffix.".matrix.gz";
print STDERR "writing matrix ($symmetricalMatrixFile)...\n" if($verbose);
writeMatrix($matrix,$matrixObject->{ inc2header },$symmetricalMatrixFile,$matrixObject->{ missingValue },$commentLine);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);