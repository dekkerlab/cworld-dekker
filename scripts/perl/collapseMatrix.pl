#!/usr/bin/perl -w

use English;
use warnings;
use strict;
use Carp qw(carp cluck croak confess);
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use POSIX qw(ceil floor);
use List::Util qw[min max];
use Cwd 'abs_path';
use Cwd;

use cworld::dekker;

my $tool=(split(/\//,abs_path($0)))[-1];

sub check_options {
    my $opts = shift;
    
    my ($inputMatrix,$verbose,$output,$collapseBy,$excludeDiagonal);
 
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
    
    if( exists($opts->{ collapseBy }) ) {
        $collapseBy = $opts->{ collapseBy };
        croak "incorrect collapseBy [$collapseBy] value (chr,name,group)" if(($collapseBy ne "chr") and ($collapseBy ne "name") and ($collapseBy ne "group"));
    } else {        
        $collapseBy="chr";
    }
    
    if( exists($opts->{ excludeDiagonal }) ) {
        $excludeDiagonal = 1;
    } else {
        $excludeDiagonal = 0;
    }
    
    $ret->{ inputMatrix }=$inputMatrix;
    $ret->{ verbose }=$verbose;
    $ret->{ output }=$output;
    $ret->{ collapseBy }=$collapseBy;
    $ret->{ excludeDiagonal }=$excludeDiagonal;
    
    return($ret,$inputMatrix,$verbose,$output,$collapseBy,$excludeDiagonal);
}

sub header2subMatrixHeader($$) {
    my $matrixObject=shift;
    my $collapseBy=shift;
    
    my $inc2header=$matrixObject->{ inc2header };
    my $header2inc=$matrixObject->{ header2inc };
    my $numYHeaders=$matrixObject->{ numYHeaders };
    my $numXHeaders=$matrixObject->{ numXHeaders };
    my $numTotalHeaders=$matrixObject->{ numTotalHeaders };
    my $missingValue=$matrixObject->{ missingValue };
    
    my %subMatrixHash=();
    
    my %special_inc2header=();
    my %special_header2inc=();
    my %special_header2subMatrixInc=();
    
    my $subMatrixIndex=0;
    for(my $i=0;$i<$numTotalHeaders;$i++) {
        my $header=$inc2header->{ xy }->{$i};
        my $subMatrix=header2subMatrix($header,$collapseBy);
        
        if(!exists($subMatrixHash{$subMatrix})) {
            $subMatrixHash{$subMatrix}=$subMatrixIndex;
            $subMatrixIndex++;
        }
    }
    
    for(my $y=0;$y<$numYHeaders;$y++) {
        my $header=$inc2header->{ y }->{$y};
        
        my $subMatrix=header2subMatrix($header,$collapseBy);
        
        croak "subMatrix [$subMatrix] does not exist!" if(!exists($subMatrixHash{$subMatrix}));
        my $subMatrixIndex=$subMatrixHash{$subMatrix};
        
        $special_inc2header{ y }{$subMatrixIndex}=$subMatrix;
        $special_header2inc{ y }{$subMatrix}=$subMatrixIndex;
        
        $special_header2subMatrixInc{ y }{$header}=$subMatrixIndex;
    }
    
    for(my $x=0;$x<$numYHeaders;$x++) {
        my $header=$inc2header->{ x }->{$x};
        
        my $subMatrix=header2subMatrix($header,$collapseBy);
        
        croak "subMatrix [$subMatrix] does not exist!" if(!exists($subMatrixHash{$subMatrix}));
        my $subMatrixIndex=$subMatrixHash{$subMatrix};
        
        $special_inc2header{ x }{$subMatrixIndex}=$subMatrix;
        $special_header2inc{ x }{$subMatrix}=$subMatrixIndex;
        
        $special_header2subMatrixInc{ x }{$header}=$subMatrixIndex;
    }
    
    return(\%special_inc2header,\%special_header2inc,\%special_header2subMatrixInc);
}

sub collapseData($$$$$;$) {
    # required
    my $inputMatrix=shift;
    my $matrixObject=shift;
    my $special_header2inc=shift;
    my $collapseBy=shift;
    my $excludeDiagonal=shift;
    # optional
    my $sigDigits=4;
    $sigDigits=shift if @_;
    
    my $header2inc=$special_header2inc;
    
    my $symmetrical=$matrixObject->{ symmetrical };
    my $headerSizing=$matrixObject->{ headerSizing };
    my $headerSpacing=$matrixObject->{ headerSpacing };
    my $missingValue=$matrixObject->{ missingValue };
    my $headerFlag=$matrixObject->{ headerFlag };
    my $verbose=$matrixObject->{ verbose };
    
    my %matrix=();
    
    my $lineNum=0;
    my @xHeaders=();
    
    print STDERR "\tgetData\n" if($verbose);
    
    my $nLines = getNumberOfLines($inputMatrix)-1;
    my $progressBucketSize=ceil($nLines / 1000);
    my $pcComplete=0;
    
    my %chrSum=();
    my %cisTrans=();
    
    my $totalInterReads=0;
    
    open(IN,inputWrapper($inputMatrix)) or croak "Could not open file [$inputMatrix] - $!";
    while(my $line = <IN>) {
        chomp($line);
        next if(($line eq "") or ($line =~ m/^#/));
        
        if($lineNum == 0) {
            @xHeaders=split(/\t/,$line);
        } else {
            my @data=split(/\t/,$line);
            my $dsize=@data;
            
            my $yHeader=$data[0];
            my $yHeadGroup=header2subMatrix($yHeader,"group");
            my $yHeadChromosome=header2subMatrix($yHeader,"liteChr");
            my $ySubMatrix=header2subMatrix($yHeader,$collapseBy);
                
            my $yIndex=-1;
            $yIndex = $header2inc->{ y }->{$yHeader} if(defined($header2inc->{ y }->{$yHeader}));
            next if($yIndex == -1);
            
            my $indexStart=1;
            my $indexEnd=$dsize;
            
            for(my $i=$indexStart;$i<$indexEnd;$i++) {
                my $cScore=$data[$i];
                
                # skip if cScore is not a valid number
                $cScore = "NA" if($cScore !~ (/^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/));
                next if($cScore eq "NA");
                
                # sparse matrix logic - do not store 0/nan (depending on quantity)
                next if( ($cScore eq $missingValue) or (($cScore ne "NA") and ($missingValue ne "NA") and ($cScore == $missingValue)) );
                
                # truncate numbers to minimal digits
                $cScore = sprintf "%.".$sigDigits."f", $cScore if($cScore ne "NA");
                
                my $xHeader=$xHeaders[$i];
                my $xHeadGroup=header2subMatrix($xHeader,"group");
                my $xHeadChromosome=header2subMatrix($xHeader,"liteChr");
                
                my $xIndex=-1;
                $xIndex = $header2inc->{ x }->{$xHeader} if(defined($header2inc->{ x }->{$xHeader}));
                next if($xIndex == -1);
                
                my $group="NA";
                $group="cis__trans" if($yHeadGroup eq $xHeadGroup);
                $group="cis__cis" if(($yHeadGroup eq $xHeadGroup) and ($yHeadChromosome eq $xHeadChromosome));
                $group="trans__trans" if($yHeadGroup ne $xHeadGroup);
                $group="trans__cis" if(($yHeadGroup ne $xHeadGroup) and ($yHeadChromosome eq $xHeadChromosome));
                
                croak "invalid group [$group] - $yHeader | $xHeader" if($group eq "NA");
                
                $group = $ySubMatrix."__".$group;
                
                my ($subMatrix,$yGroup,$xGroup)=split(/__/,$group);
                $cisTrans{$subMatrix}{$yGroup}{$xGroup} += $cScore;
                
                next if(($excludeDiagonal) and ($yHeadChromosome eq $xHeadChromosome));
                next if(($excludeDiagonal) and ($yIndex == $xIndex));
                
                $matrix{$yIndex}{$xIndex} += $cScore;
                
                $chrSum{$yIndex} += $cScore if($yIndex != $xIndex);
                $totalInterReads += $cScore if($yIndex != $xIndex);
                
            }
        }
                
        $pcComplete = 100 if($lineNum == ($nLines-1));
        print STDERR "\e[A" if(($verbose) and ($lineNum != 0));
        printf STDERR "\t%.2f%% complete ($lineNum/".($nLines-1).")...\n", $pcComplete if($verbose);
        $pcComplete = round((($lineNum/$nLines)*100),2);
        $lineNum++;
    }
    close(IN);
    
    $pcComplete=100;
    print STDERR "\e[A" if($verbose);
    printf STDERR "\t%.2f%% complete ($lineNum/$nLines)...\n", $pcComplete if($verbose);

    return(\%matrix,\%cisTrans,\%chrSum,$totalInterReads);
}



sub intro() {
    print STDERR "\n";
    
    print STDERR "Tool:\t\t".$tool."\n";
    print STDERR "Version:\t".$cworld::dekker::VERSION."\n";
    print STDERR "Summary:\tcollapse matrix by (chr,name,group), sum signal\n";
    
    print STDERR "\n";
}

sub help() {
    intro();
    
    print STDERR "Usage: perl collapseMatrix.pl [OPTIONS] -i <inputMatrix>\n";
    
    print STDERR "\n";
    
    print STDERR "Required:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-i", "[]", "input matrix file");
    
    print STDERR "\n";
    
    print STDERR "Options:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-v", "[]", "FLAG, verbose mode");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-o", "[]", "prefix for output file(s)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--cb", "[chr]", "method to collapse headers (chr,name,group)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ed", "[]", "FLAG, exclude diagonal");
    
    print STDERR "\n";
    
    print STDERR "Notes:";
    print STDERR "
    This script collapses headers by chr,name,group.
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

my %options;
my $results = GetOptions( \%options,'inputMatrix|i=s','verbose|v','output|o=s','collapseBy|cb=s','excludeDiagonal|ed') or croak help();
my ($ret,$inputMatrix,$verbose,$output,$collapseBy,$excludeDiagonal)=check_options( \%options );

intro() if($verbose);

#get the absolute path of the script being used
my $cwd = getcwd();
my $fullScriptPath=abs_path($0);
my @fullScriptPathArr=split(/\//,$fullScriptPath);
@fullScriptPathArr=@fullScriptPathArr[0..@fullScriptPathArr-3];
my $scriptPath=join("/",@fullScriptPathArr);

croak "inputMatrix [$inputMatrix] does not exist" if(!(-e $inputMatrix));

# get matrix information
my $matrixObject=getMatrixObject($inputMatrix,$output,$verbose);
my $inc2header=$matrixObject->{ inc2header };
my $header2inc=$matrixObject->{ header2inc };
my $numYHeaders=$matrixObject->{ numYHeaders };
my $numXHeaders=$matrixObject->{ numXHeaders };
my $missingValue=$matrixObject->{ missingValue };
my $symmetrical=$matrixObject->{ symmetrical };
my $inputMatrixName=$matrixObject->{ inputMatrixName };
$output=$matrixObject->{ output };

print STDERR "parsing headers ..\n" if($verbose);
my ($special_inc2header,$special_header2inc,$special_header2subMatrixInc)=header2subMatrixHeader($matrixObject,$collapseBy);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

#read Matrix
print STDERR "collapsing data ...\n" if($verbose);
my ($collapsedMatrix,$cisTrans,$chrSum,$totalInterReads)=collapseData($inputMatrix,$matrixObject,$special_header2subMatrixInc,$collapseBy,$excludeDiagonal);
print STDERR "\tdone.\n" if($verbose);

$matrixObject->{ header2inc }=$special_header2inc;
$matrixObject->{ inc2header }=$special_inc2header;

print STDERR "\n" if($verbose);

print STDERR "writing raw collapsed matrix ...\n" if($verbose);
my $collapsedMatrixFile=$output.".".$collapseBy.".collapsed.matrix.gz";
writeMatrix($collapsedMatrix,$special_inc2header,$collapsedMatrixFile,"NA");
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

my %possibleGroups = (
            cis => 1,
            trans => 1,
    );

my $cisTransFile=$output.".".$collapseBy.".cisTrans.txt";
open(OUT,outputWrapper($cisTransFile)) or croak "Could not open file [$cisTransFile] - $!";
print OUT "subMatrix\tcis__cis\tcis__trans\ttrans__cis\ttrans__trans\n";

foreach my $subMatrix ( sort keys %{$cisTrans} ) {
    print OUT "$subMatrix";
    foreach my $yGroup ( sort keys %possibleGroups ) {
        foreach my $xGroup ( sort keys %possibleGroups ) {
            my $sum=0;
            $sum=$cisTrans->{$subMatrix}->{$yGroup}->{$xGroup} if(exists($cisTrans->{$subMatrix}->{$yGroup}->{$xGroup}));
            print OUT "\t$sum";
        }
    }
    print OUT "\n";
}
close(OUT);

print STDERR "plotting cisTrans data ...\n" if($verbose);
system("Rscript '".$scriptPath."/R/cisTrans.R' '".$cisTransFile."' > /dev/null");
print STDERR "\n" if($verbose);

print STDERR "totalInterReads = $totalInterReads\n" if($verbose);

print STDERR "\n" if($verbose);

my $normalizedMatrix={};
print STDERR "normalizing matrix ...\n" if($verbose);
my $sumPCReads=0;

my $numYSubMatrices = keys %{$special_inc2header->{ y }};
my $numXSubMatrices = keys %{$special_inc2header->{ x }};

for(my $ySubIndex=0;$ySubIndex<$numYSubMatrices;$ySubIndex++) {
    my $ySubHeader=$special_inc2header->{ y }->{$ySubIndex};
    my $yInterSum=$chrSum->{$ySubIndex};
    
    for(my $xSubIndex=0;$xSubIndex<$numXSubMatrices;$xSubIndex++) {
        my $xSubHeader=$special_inc2header->{ x }->{$xSubIndex};
        my $xInterSum=$chrSum->{$xSubIndex};
        
        my $score="NA";
        $score=$collapsedMatrix->{$ySubIndex}->{$xSubIndex} if(exists($collapsedMatrix->{$ySubIndex}->{$xSubIndex}));
        $score = "NA" if(($score ne "NA") and ($score eq "NA"));
        
        my $logRatio="NA";
        if($score ne "NA") {
            # rachel expected calculation
            my $numerator = $score;
            my $denominator = 0;
            $denominator = ( ( (($yInterSum/$totalInterReads) * ($xInterSum/($totalInterReads-$yInterSum))) + (($xInterSum/$totalInterReads)*($yInterSum/($totalInterReads-$xInterSum))) ) * ($totalInterReads/2) ) if(($totalInterReads != 0) and (($totalInterReads-$yInterSum) != 0));
        
            my $normalizedScore = 0;
            $normalizedScore = ($numerator / $denominator) if($denominator != 0);
            
            $logRatio = (log($normalizedScore)/log(2)) if($normalizedScore != 0) and (($normalizedScore) != 0);
            
            $logRatio = "NA" if($ySubHeader eq $xSubHeader);
        }
        
        $normalizedMatrix->{$ySubIndex}->{$xSubIndex} = $logRatio;
        
    }
}
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

print STDERR "writing normalized collapsed matrix ...\n" if($verbose);
my $normalizedCollapsedMatrixFile=$output.".".$collapseBy.".collapsed.normalized.matrix.gz";
writeMatrix($normalizedMatrix,$special_inc2header,$normalizedCollapsedMatrixFile,0);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);