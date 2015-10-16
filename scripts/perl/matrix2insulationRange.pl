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

    my ($inputMatrix,$verbose,$output,$insulationSquareStart,$insulationSquareEnd,$insulationSquareStep,$insulationMode,$minDistance,$maxDistance,$excludeZero);
    
    my $ret={};
    
    if( exists($opts->{ inputMatrix }) ) {
        $inputMatrix = $opts->{ inputMatrix };
    } else {
        print STDERR "\nERROR: Option inputMatrix|i is required.\n";
        help();
    }
    
    if( exists($opts->{ verbose }) ) {
        $verbose = 1;
    } else {
        $verbose = 0;
    }
    
    if( exists($opts->{ output }) ) {
        $output = $opts->{ output };
    } else {
        $output = "";
    }
    
    if( exists($opts->{ insulationSquareStart }) ) {
        $insulationSquareStart = $opts->{ insulationSquareStart };
    } else {
        $insulationSquareStart=200000;
    }
    
    if( exists($opts->{ insulationSquareEnd }) ) {
        $insulationSquareEnd = $opts->{ insulationSquareEnd };
    } else {
        $insulationSquareEnd=500000;
    }
    
    if( exists($opts->{ insulationSquareStep }) ) {
        $insulationSquareStep = $opts->{ insulationSquareStep };
    } else {
        $insulationSquareStep=10000;
    }
    
    if( exists($opts->{ insulationMode }) ) {
        $insulationMode = $opts->{ insulationMode };
    } else {
        $insulationMode="mean";
    }
    
    if( exists($opts->{ minDistance }) ) {
        $minDistance = $opts->{ minDistance };
    } else {
        $minDistance=undef;
    }
    
    if( exists($opts->{ maxDistance }) ) {
        $maxDistance = $opts->{ maxDistance };
    } else {
        $maxDistance = undef;
    }
    
    if( exists($opts->{ excludeZero }) ) {
        $excludeZero = 1;
    } else {
        $excludeZero = 0;
    }
    
    $ret->{ inputMatrix }=$inputMatrix;
    $ret->{ verbose }=$verbose;
    $ret->{ output }=$output;
    $ret->{ insulationSquareStart }=$insulationSquareStart;
    $ret->{ insulationSquareEnd }=$insulationSquareEnd;
    $ret->{ insulationSquareStep }=$insulationSquareStep;
    $ret->{ insulationMode }=$insulationMode;
    $ret->{ minDistance }=$minDistance;
    $ret->{ maxDistance }=$maxDistance;
    $ret->{ excludeZero }=$excludeZero;
    
    return($ret,$inputMatrix,$verbose,$output,$insulationSquareStart,$insulationSquareEnd,$insulationSquareStep,$insulationMode,$minDistance,$maxDistance,$excludeZero);
}

sub intro() {
    print STDERR "\n";
    
    print STDERR "Tool:\t\t".$tool."\n";
    print STDERR "Version:\t".$cworld::dekker::VERSION."\n";
    print STDERR "Summary:\tcalculate insulation index (TADs) over range of square sizes\n";
    
    print STDERR "\n";
}

sub help() {
    intro();
    
    print STDERR "Usage: perl matrix2insulationRange.pl [OPTIONS] -i <inputMatrix>\n";
    
    print STDERR "\n";
    
    print STDERR "Required:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-i", "[]", "input matrix file");
    
    print STDERR "\n";
    
    print STDERR "Options:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-v", "[]", "FLAG, verbose mode");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-o", "[]", "prefix for output file(s)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--im", "[]", "optional, insulation mode (how to aggregrate signal within insulation square), mean,sum,median");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--istart", "[]", "insulation square size start, minumum bin size");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--iend", "[]", "insulation square size end, maximum bin size");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--istep", "[]", "insulation square size step, step size (in bp) for insulation square range");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--minDist", "[]", "minimum allowed interaction distance, exclude < N distance (in BP)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--maxDist", "[]", "maximum allowed interaction distance, exclude > N distance (in BP)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ez", "[]", "FLAG, ignore 0s in all calculations");    

    print STDERR "\n";
    
    print STDERR "Notes:";
    print STDERR "
    This script calculates the insulation index over a range of insulation square sizes.
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


sub calculateInsulation($$$$$$$$) {
    my $matrixObject=shift;
    my $matrix=shift;
    my $insulationSquareStart=shift;
    my $insulationSquareEnd=shift;
    my $insulationSquareStep=shift;
    my $insulationMode=shift;
    my $excludeZero=shift;
    my $scriptPath=shift;
    
    my $inc2header=$matrixObject->{ inc2header };
    my $header2inc=$matrixObject->{ header2inc };
    my $numHeaders=$matrixObject->{ numYHeaders };
    my $headerSizing=$matrixObject->{ headerSizing };
    my $headerSpacing=$matrixObject->{ headerSpacing };
    my $verbose=$matrixObject->{ verbose };
    
    my %matrixInsulation=();
    
    my $pcComplete=0;
    
    for(my $y=0;$y<$numHeaders;$y++) {
        my $yHead=$inc2header->{ y }->{$y};
        
        my $tmp_squareEndBins=ceil(($insulationSquareEnd-($headerSizing-$headerSpacing))/$headerSpacing)+1;
        
        my @boxData=();
        my $count=0;
        my $sum=0;
        my $skipFlag=0;
        my $nullCount=0;
        
        # for all insulation square sizes
        for(my $offset=1;$offset<=$tmp_squareEndBins;$offset++) {
            # offset = num bins to reach
            
            my $expectedDataCount=($offset*$offset);
            
            # col
            my $col_yEnd=($y-$offset);
            last if($col_yEnd < 0);
            my $col_yStart=($col_yEnd+$offset)-1;
            last if($col_yStart > $numHeaders);
            my $col_x=($y+$offset);
            last if($col_x > ($numHeaders-1));
            my $tmpColumnHeader=$inc2header->{ x }->{$col_x};
            
            # row 
            my $row_xStart=($y+1);
            my $row_xEnd=($row_xStart+$offset)-1;
            last if($row_xEnd > $numHeaders);
            my $row_y=($y-$offset);
            last if($row_y < 0);
            my $tmpRowHeader=$inc2header->{ y }->{$row_y};
            
            croak "col/row out of alignment!" if(($col_yEnd != $row_y) or ($col_x != $row_xEnd));
            
            my $corner_y=$row_y if($col_yEnd == $row_y);
            my $corner_x=$col_x  if($col_x == $row_xEnd);
            
            # get the column
            my  @tmpColumn=();
            for(my $tmp_y=$col_yStart;$tmp_y>$col_yEnd;$tmp_y--) {
                
                my $inten=$matrixObject->{ missingValue };
                $inten=$matrix->{$tmp_y}->{$col_x} if(exists($matrix->{$tmp_y}->{$col_x}));
                $skipFlag = 1 if($inten eq "NA");
                $nullCount++ if($inten eq "NA");
                next if($inten eq "NA");
                next if(($inten eq 0) and ($excludeZero));
                    
                push(@boxData,$inten);
                push(@tmpColumn,$inten);
                $sum+=$inten;
                $count++;
            }
            my $tmpColumnStats=listStats(\@tmpColumn) if(@tmpColumn > 0);
            my $tmpColumnAggregate="NA";
            $tmpColumnAggregate=$tmpColumnStats->{ iqrMean } if(defined($tmpColumnStats->{ iqrMean }));
            
            # only do row if offset > 1 (to avoid double counting 1x1 square)
            if($offset > 1) {
                # get the row
                my @tmpRow=();
                for(my $tmp_x=$row_xStart;$tmp_x<$row_xEnd;$tmp_x++) {
                    
                    my $inten=$matrixObject->{ missingValue };
                    $inten=$matrix->{$row_y}->{$tmp_x} if(exists($matrix->{$row_y}->{$tmp_x}));
                    $skipFlag = 1 if($inten eq "NA");
                    $nullCount++ if($inten eq "NA");
                    next if($inten eq "NA");
                    next if(($inten eq 0) and ($excludeZero));
                        
                    push(@boxData,$inten);
                    push(@tmpRow,$inten);
                    $sum+=$inten;
                    $count++;
                }
                my $tmpRowStats=listStats(\@tmpRow) if(@tmpRow > 0);
                my $tmpRowAggregate="NA";
                $tmpRowAggregate=$tmpRowStats->{ iqrMean } if(defined($tmpRowStats->{ iqrMean }));
            }
            
            # get the corner
            my $inten=$matrixObject->{ missingValue };
            $inten=$matrix->{$corner_y}->{$corner_x} if(exists($matrix->{$corner_y}->{$corner_x}));
            $skipFlag = 1 if($inten eq "NA");
            $nullCount++ if($inten eq "NA");
            if(!(($inten eq "NA") or (($inten eq 0) and ($excludeZero)))) {
                push(@boxData,$inten);
                $sum+=$inten;
                $count++;
            } 
            
            # dump current aggregrate of insulation square
            
            my $boxDataSize=scalar @boxData;            
            next if($boxDataSize == 0);
            next if(($skipFlag) and ($insulationMode eq "sum"));
            
            my $mean="NA";
            $mean = ($sum/$count) if($count != 0);
            
            my $boxDataAggregate="NA";
            $boxDataAggregate=$sum if($insulationMode eq "sum");
            $boxDataAggregate=$mean if($insulationMode eq "mean");
                        
            my $boxDataStats=listStats(\@boxData) if(($insulationMode ne "sum") and ($insulationMode ne "mean"));
            $boxDataAggregate=$boxDataStats->{ $insulationMode } if(($insulationMode ne "sum") and ($insulationMode ne "mean"));
            my $boxDataZeroPercent=$boxDataStats->{ zeroPC } if(($insulationMode ne "sum") and ($insulationMode ne "mean"));
                        
            my $tmp_square=($offset * $headerSpacing)+($headerSizing-$headerSpacing);
            
            # if > half values are NA, refuse to calculate insulation value
            next if($nullCount > ($expectedDataCount/2));
            
            # store result in hash header->insulation
            $matrixInsulation{$tmp_square}{$yHead}=$boxDataAggregate;
            
        }
        
        $pcComplete = 100 if($y == ($numHeaders-1));
        print STDERR "\e[A" if(($verbose) and ($y != 0));
        printf STDERR "\t%.2f%% complete ($y/".($numHeaders-1).")...\n", $pcComplete if($verbose);
        $pcComplete = round((($y/($numHeaders-1))*100),2);
            
    }
    
    return(\%matrixInsulation);
    
}

sub normalizeInsulation($$) {
    my $matrixObject=shift;
    my $matrixInsulation=shift;
    
    my $inc2header=$matrixObject->{ inc2header };
    my $header2inc=$matrixObject->{ header2inc };
    my $numHeaders=$matrixObject->{ numYHeaders };
    my $headerSizing=$matrixObject->{ headerSizing };
    my $headerSpacing=$matrixObject->{ headerSpacing };
    my $verbose=$matrixObject->{ verbose };
    
    my %matrixNormalizedInsulation=();
    
    foreach my $square (sort { $a <=> $b } keys %$matrixInsulation) {
        my @insulationSignals=();
        for(my $y=0;$y<$numHeaders;$y++) {
            my $yHead=$inc2header->{ y }->{$y};
            
            next if(!exists($matrixInsulation->{$square}->{$yHead}));
            
            my $insulation="NA";
            $insulation=$matrixInsulation->{$square}->{$yHead} if(exists($matrixInsulation->{$square}->{$yHead}));
            next if($insulation eq "NA");
            
            push(@insulationSignals,$insulation)
        }
    
        croak "no available data points!" if(@insulationSignals == 0);
    
        my $tmpStats=listStats(\@insulationSignals);
        my $meanInsulation=$tmpStats->{ mean };
        print STDERR "\t$square\tmean centering data ($meanInsulation)\n" if($verbose);
    
        for(my $y=0;$y<$numHeaders;$y++) {
            my $yHead=$inc2header->{ y }->{$y};
        
            next if(!exists($matrixInsulation->{$square}->{$yHead}));
            
            my $insulation="NA";
            $insulation=$matrixInsulation->{$square}->{$yHead} if(exists($matrixInsulation->{$square}->{$yHead}));
            
            $matrixNormalizedInsulation{$square}{$yHead}=$insulation;
            
            #mean center the insulation scores
            my $normalizedInsulation = "NA";
            $normalizedInsulation = (log($insulation/$meanInsulation)/log(2)) if(($meanInsulation ne "NA") and ($insulation ne "NA") and ($meanInsulation != 0) and ($insulation != 0));
            
            $matrixNormalizedInsulation{$square}{$yHead}=$normalizedInsulation;
        }
    }
    
    return(\%matrixNormalizedInsulation);
}


sub outputInsulation($$$$) {
    my $matrixObject=shift;
    my $matrixInsulation=shift;
    my $insulationRangeFileName=shift;
    my $commentLine=shift;
    
    my $inc2header=$matrixObject->{ inc2header };
    my $header2inc=$matrixObject->{ header2inc };
    my $numHeaders=$matrixObject->{ numYHeaders };
    my $headerSizing=$matrixObject->{ headerSizing };
    my $headerSpacing=$matrixObject->{ headerSpacing };
    
    # main insulation file
    open(OUT,outputWrapper($insulationRangeFileName,$commentLine)) or croak "Could not open file [$insulationRangeFileName] - $!";
    
    for(my $y=0;$y<$numHeaders;$y++) {
        my $yHead=$inc2header->{ y }->{$y};
        my $yHeadObject=getHeaderObject($yHead);
        my $yHeadChromosome=$yHeadObject->{ chromosome };
        my $yHeadStart=$yHeadObject->{ start };
        my $yHeadEnd=$yHeadObject->{ end };
        my $yHeadMidpoint=round(($yHeadStart+$yHeadEnd)/2);
        
        print OUT "\tmidpoint-$yHeadMidpoint";
    }
    
    print OUT "\n";
            
    foreach my $square (sort { $b <=> $a } keys %$matrixInsulation) {
        print OUT "square-$square";
        
        for(my $y=0;$y<$numHeaders;$y++) {
            my $yHead=$inc2header->{ y }->{$y};
            
            my $insulation="NA";
            $insulation=$matrixInsulation->{$square}->{$yHead} if(exists($matrixInsulation->{$square}->{$yHead}));
        
            print OUT "\t$insulation";
        }
        print OUT "\n";
    }
    
    close(OUT);
    
    return($insulationRangeFileName);
}

my %options;
my $results = GetOptions( \%options,'inputMatrix|i=s','verbose|v','output|o=s','insulationSquareStart|istart=i','insulationSquareEnd|iend=i','insulationSquareStep|istep=i','insulationMode|im=s','minDistance|minDist=i','maxDistance|maxDist=i','excludeZero|ez') or croak help();
my ($ret,$inputMatrix,$verbose,$output,$insulationSquareStart,$insulationSquareEnd,$insulationSquareStep,$insulationMode,$minDistance,$maxDistance,$excludeZero)=check_options( \%options );

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
my $numHeaders=$matrixObject->{ numTotalHeaders };
my $missingValue=$matrixObject->{ missingValue };
my $contigList=$matrixObject->{ contigList };
my $numContigs=$matrixObject->{ numContigs };
my $index2contig=$matrixObject->{ index2contig };
my $symmetrical=$matrixObject->{ symmetrical };
my $equalHeaderFlag=$matrixObject->{ equalHeaderFlag };
my $inputMatrixName=$matrixObject->{ inputMatrixName };
$output=$matrixObject->{ output };

croak "matrix must be symmetrical" if($symmetrical == 0);

my $imageWidth=$numHeaders*2;
$imageWidth=900 if($imageWidth < 900);

print STDERR "contigs:\n" if($verbose);
for(my $c=0;$c<$numContigs;$c++) {
    my $contig=$index2contig->{ xy }->{$c};
    my $contigStart=$contigList->{$contig}->{ contigStart };
    my $contigEnd=$contigList->{$contig}->{ contigEnd };
    my $contigAssembly=$contigList->{$contig}->{ contigAssembly };
    my $contigChromosome=$contigList->{$contig}->{ contigChromosome };
    my $contigLength=$contigList->{$contig}->{ contigLength };
    print STDERR "\t$contig\t$contigAssembly\t$contigChromosome\t$contigStart - $contigEnd [$contigLength]\n" if($verbose);
}
    
# get fragment spacing (i.e. bin size)
my ($equalSpacingFlag,$equalSizingFlag,$headerSpacing,$headerSizing)=getHeaderSpacing($inc2header->{ y });

carp "insulationSquareStart ($insulationSquareStart) < headerSpacing ($headerSpacing), setting to [$headerSpacing]" if($insulationSquareStart < $headerSpacing);
$insulationSquareStart=$headerSpacing if($insulationSquareStart < $headerSpacing);
carp "insulationSquareEnd ($insulationSquareEnd) < headerSpacing ($headerSpacing), setting to [$headerSpacing]" if($insulationSquareEnd < $headerSpacing);
$insulationSquareEnd=$headerSpacing if($insulationSquareEnd < $headerSpacing);
carp "insulationSquareStep ($insulationSquareStep) < headerSpacing ($headerSpacing), setting to [$headerSpacing]" if($insulationSquareStep < $headerSpacing);
$insulationSquareStep=$headerSpacing if($insulationSquareStep < $headerSpacing);

print STDERR "\n" if($verbose);

# add suffix
$output .= "--istart".$insulationSquareStart."--iend".$insulationSquareEnd."--istep".$insulationSquareStep."--im".$insulationMode;

carp "insulation image size is < numHeaders ($numHeaders), visual artifacts are possible!" if($numHeaders > $imageWidth);

#read Matrix
my ($matrix)=getData($inputMatrix,$matrixObject,$verbose,(($insulationSquareEnd)+$headerSizing));

# reset sparse value to NA
$matrixObject->{ missingValue }="NA";
$missingValue=$matrixObject->{ missingValue };

print STDERR "\n" if($verbose);

my $subsetMatrixFile=$output.".subset.matrix.gz";
print STDERR "Writing matrix to file ($subsetMatrixFile)...\n" if($verbose);
writeMatrix($matrix,$inc2header,$subsetMatrixFile,$missingValue);
print STDERR "\tcomplete\n" if($verbose);

print STDERR "\n" if($verbose);

# calculate the insulation index for each bin and store in a new data struct.
print STDERR "calculating insulation index [$insulationSquareStart - $insulationSquareEnd] [$insulationSquareStep] ...\n" if($verbose);
my ($matrixInsulation)=calculateInsulation($matrixObject,$matrix,$insulationSquareStart,$insulationSquareEnd,$insulationSquareStep,$insulationMode,$excludeZero,$scriptPath);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

# mean-center insulation data
print STDERR "mean-centering insulation index...\n" if($verbose);
my ($normalizedInsulation)=normalizeInsulation($matrixObject,$matrixInsulation);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

# write insulation data to file
my $insulationRangeFileName=$output.".insulation.matrix.gz";
print STDERR "outputing insulation...\n" if($verbose);
($insulationRangeFileName)=outputInsulation($matrixObject,$normalizedInsulation,$insulationRangeFileName,$commentLine);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);
