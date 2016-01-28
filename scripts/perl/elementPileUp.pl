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

    my ($inputMatrix,$elementBedFile,$verbose,$output,$loessFile,$elementZoneSize,$aggregrateMode,$minDistance,$maxDistance,$excludeZero);
    
    my $ret={};
    
    if( exists($opts->{ inputMatrix }) ) {
        $inputMatrix = $opts->{ inputMatrix };
    } else {
        print STDERR "\nERROR: Option inputMatrix|i is required.\n";
        help();
    }
    
    if( exists($opts->{ elementBedFile }) ) {
        $elementBedFile = $opts->{ elementBedFile };
    } else {
        print STDERR "\nERROR: Option elementBedFile|ebf is required.\n";
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
    
    if( exists($opts->{ elementZoneSize }) ) {
        $elementZoneSize = $opts->{ elementZoneSize };
    } else {
        $elementZoneSize=100000;
    }
    
    if( exists($opts->{ aggregrateMode }) ) {
        $aggregrateMode = $opts->{ aggregrateMode };
    } else {
        $aggregrateMode="mean";
    }
    
    if( exists($opts->{ minDistance }) ) {
        $minDistance = $opts->{ minDistance };
    } else {
        $minDistance=undef;
    }
    
    if( exists($opts->{ maxDistance }) ) {
        $maxDistance = $opts->{ maxDistance };
    } else {
        $maxDistance=100000;
    }
    
    if( exists($opts->{ excludeZero }) ) {
        $excludeZero = 1;
    } else {
        $excludeZero = 0;
    }
    
    $ret->{ inputMatrix }=$inputMatrix;
    $ret->{ elementBedFile }=$elementBedFile;
    $ret->{ verbose }=$verbose;
    $ret->{ output }=$output;
    $ret->{ loessFile }=$loessFile;
    $ret->{ elementZoneSize }=$elementZoneSize;
    $ret->{ aggregrateMode }=$aggregrateMode;
    $ret->{ minDistance }=$minDistance;
    $ret->{ maxDistance }=$maxDistance;
    $ret->{ excludeZero }=$excludeZero;
    
    return($ret,$inputMatrix,$elementBedFile,$verbose,$output,$loessFile,$elementZoneSize,$aggregrateMode,$minDistance,$maxDistance,$excludeZero);
}

sub intro() {
    print STDERR "\n";
    
    print STDERR "Tool:\t\t".$tool."\n";
    print STDERR "Version:\t".$cworld::dekker::VERSION."\n";
    print STDERR "Summary:\tpile up cData around specified list of 'elements'\n";
    
    print STDERR "\n";
}

sub help() {
    intro();
    
    print STDERR "Usage: perl elementPileUp.pl [OPTIONS] -i <inputMatrix> -ebf <elementBedFile>\n";
    
    print STDERR "\n";
    
    print STDERR "Required:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-i", "[]", "input matrix file");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ebf", "[]", "element bed file");
    
    print STDERR "\n";
    
    print STDERR "Options:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-v", "[]", "FLAG, verbose mode");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-o", "[]", "prefix for output file(s)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ezs", "[]", "elementZoneSize, size of zone to use around element (x axis - in BP)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--am", "[]", "aggregrate mode, how to aggregrate signal [mean,sum,median]");    
    printf STDERR ("\t%-10s %-10s %-10s\n", "--minDist", "[]", "minimum allowed interaction distance, exclude < N distance (in BP)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--maxDist", "[]", "maximum allowed interaction distance, exclude > N distance (in BP)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ez", "[]", "FLAG, ignore 0s in all calculations");    
    

    print STDERR "\n";
    
    print STDERR "Notes:";
    print STDERR "
    This script aggregrates cData surrounding a list of 'elements'.
    Input Matrix can be TXT or gzipped TXT.
    ElementListFile must be BED5+ format.
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

sub elementPileUp($$$$$$$$$$) {
    my $matrixObject=shift;
    my $matrix=shift;
    my $elementFileName=shift;
    my $elements=shift;
    my $elementZoneSize=shift;
    my $minDistance=shift;
    my $maxDistance=shift;
    my $aggregrateMode=shift;
    my $excludeZero=shift;
    my $commentLine=shift;
    
    my $inc2header=$matrixObject->{ inc2header };
    my $header2inc=$matrixObject->{ header2inc };
    my $numXHeaders=$matrixObject->{ numXHeaders };
    my $numYHeaders=$matrixObject->{ numYHeaders };
    my $numTotalHeaders=$matrixObject->{ numTotalHeaders };
    my $headerSizing=$matrixObject->{ headerSizing };
    my $headerSpacing=$matrixObject->{ headerSpacing };
    my $output=$matrixObject->{ output };
    my $verbose=$matrixObject->{ verbose };
    
    my $maxDistanceLimit=($numTotalHeaders * $headerSpacing)+($headerSizing-$headerSpacing);
    my $maxDistanceLimit_bins=$numTotalHeaders;
    
    # if distance limit = NA, set to max
    $maxDistance=$maxDistanceLimit if($maxDistance eq "NA");
    
    # do not allow either element/distance to reach beyond matrix bounds
    $elementZoneSize=$maxDistanceLimit if($elementZoneSize > $maxDistanceLimit);
    $maxDistance=$maxDistanceLimit if($maxDistance > $maxDistanceLimit);
    
    my $elementZoneSize_bins=ceil(($elementZoneSize-($headerSizing-$headerSpacing))/$headerSpacing);
    my $elementZoneSize_binsDistance=($elementZoneSize_bins * $headerSpacing)+($headerSizing-$headerSpacing);
    print STDERR "\telementZoneSize\t$elementZoneSize\t$elementZoneSize_bins\t$elementZoneSize_binsDistance\n" if($verbose);

    my $maxDistance_bins=ceil(($maxDistance-($headerSizing-$headerSpacing))/$headerSpacing);
    my $maxDistance_binsDistance=($maxDistance_bins * $headerSpacing)+($headerSizing-$headerSpacing);
    print STDERR "\tmaxDistance\t$maxDistance\t$maxDistance_bins\t$maxDistance_binsDistance\n" if($verbose);
    
    my $pileUpMatrix={};
    
    my $numElements=@{$elements};
    my $pcComplete=0;
    for(my $a=0;$a<$numElements;$a++) {
        my $element=$elements->[$a]->{ name };
        my $elementObject=getHeaderObject($element,1);
        
        my $elementIndex=-1;
        $elementIndex=$header2inc->{ xy }->{$element} if(exists($header2inc->{ xy }->{$element}));
        
        croak "invalid header! [$element]" if($elementIndex == -1);
        
        my $leftBoundIndex=$elementIndex-$elementZoneSize_bins;
        my $rightBoundIndex=$elementIndex+$elementZoneSize_bins;
        
        for(my $x=-$elementZoneSize_bins;$x<=$elementZoneSize_bins;$x++) {
            my $origX=($elementIndex+$x);
            my $tmpHeader=$inc2header->{ xy }->{$origX};
            for(my $y=1;$y<=$maxDistance_bins;$y++) {
                
                my $tmpX=($origX+$y);
                my $tmpY=($origX-$y);
                
                # catch out of bounds
                next if($tmpX >= $numTotalHeaders);
                next if($tmpX < 0);
                next if($tmpY >= $numTotalHeaders);
                next if($tmpY < 0);
                
                my $tmp_xHeader=$inc2header->{ xy }->{$tmpX};
                my $tmp_yHeader=$inc2header->{ xy }->{$tmpY};
                
                my $tmp_xHeaderObject=getHeaderObject($tmp_xHeader,1);
                my $tmp_yHeaderObject=getHeaderObject($tmp_yHeader,1);
                my $interactionDistance=getInteractionDistance($matrixObject,$tmp_yHeaderObject,$tmp_xHeaderObject,1);
                
                next if($interactionDistance == -1);
                
                my $inten=$matrixObject->{ missingValue };
                $inten=$matrix->{$tmpY}->{$tmpX} if(exists($matrix->{$tmpY}->{$tmpX}));
                $inten = "NA" if($inten eq "NA");
                $inten = "NA" if(($tmpX > $numTotalHeaders) || ($tmpY > $numTotalHeaders));

                next if($inten eq "NA");
                
                if(($aggregrateMode ne "sum") and ($aggregrateMode ne "mean")) {
                    @{$pileUpMatrix->{$y}->{$x}}=() if(!exists($pileUpMatrix->{$y}->{$x}));
                    push(@{$pileUpMatrix->{$y}->{$x}},$inten);
                } else {
                    $pileUpMatrix->{$y}->{$x}->{ sum }+=$inten;
                    $pileUpMatrix->{$y}->{$x}->{ count }++;
                }
            }
        }
        
        $pcComplete = 100 if($a == ($numElements-1));
        print STDERR "\e[A" if(($verbose) and ($a != 0));
        printf STDERR "\t%.2f%% complete ($a/".($numElements-1).")...\n", $pcComplete if($verbose);
        $pcComplete = round((($a/($numElements-1))*100),2);
        
    }
    
    print STDERR "\n" if($verbose);
    
    my $pileUpMatrixFile=$output."___".$elementFileName.".pileUp.matrix.gz";
    print STDERR "writing pileUpMatrix ...\n" if($verbose);
    
    open(OUT,outputWrapper($pileUpMatrixFile,$commentLine)) or croak "Could not open file [$pileUpMatrixFile] - $!";
    
    for(my $x=-$elementZoneSize_bins;$x<=$elementZoneSize_bins;$x++) {
        my $xLabel=($x*$headerSpacing)."bp";
        $xLabel="+".$xLabel if($x > 0);
        print OUT "\t$xLabel";
    }
    
    print OUT "\n";
    
    for(my $y=$maxDistance_bins;$y>0;$y--) {
        my $yLabel="+".($y*$headerSpacing)."bp";
        print OUT "$yLabel";
        for(my $x=-$elementZoneSize_bins;$x<=$elementZoneSize_bins;$x++) {
            
            my $score="NA";
            if(($aggregrateMode ne "sum") and ($aggregrateMode ne "mean")) {
                my @tmpArr=();
                @tmpArr=@{$pileUpMatrix->{$y}->{$x}} if(exists($pileUpMatrix->{$y}->{$x}));
                my $tmpArrSize=scalar @tmpArr;
                my $tmpArrStats=listStats(\@tmpArr) if(@tmpArr > 0);
                $score=$tmpArrStats->{ $aggregrateMode } if(exists($tmpArrStats->{ $aggregrateMode }));
            } else {
                my $sum="NA";
                $sum=$pileUpMatrix->{$y}->{$x}->{ sum } if(exists($pileUpMatrix->{$y}->{$x}->{ sum }));
                my $count="NA";
                $count=$pileUpMatrix->{$y}->{$x}->{ count } if(exists($pileUpMatrix->{$y}->{$x}->{ count }));
                my $mean="NA";
                $mean=($sum/$count) if((($count ne "NA") and ($sum ne "NA")) and ($count != 0));
                $score=$sum if($aggregrateMode eq "sum");
                $score=$mean if($aggregrateMode eq "mean");
            }
            print OUT "\t$score";
        }
        print OUT "\n";
    }
    
    close(OUT);
    
    print STDERR "\tdone\n" if($verbose);
    
    return($pileUpMatrixFile);
    
}

my %options;
my $results = GetOptions( \%options,'inputMatrix|i=s','elementBedFile|ebf=s','verbose|v','output|o=s','loessObjectFile|lof=s','elementZoneSize|ezs=i','aggregrateMode|am=s','minDistance|minDist=i','maxDistance|maxDist=i','excludeZero|ez') or croak help();
my ($ret,$inputMatrix,$elementBedFile,$verbose,$output,$loessObjectFile,$elementZoneSize,$aggregrateMode,$minDistance,$maxDistance,$excludeZero)=check_options( \%options );

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
my $missingValue=$matrixObject->{ missingValue };
my $symmetrical=$matrixObject->{ symmetrical };
my $inputMatrixName=$matrixObject->{ inputMatrixName };
$output=$matrixObject->{ output };

croak "matrix must be symmetrical!" if($symmetrical == 0);

print STDERR "validating [$elementBedFile] ...\n" if($verbose);
validateBED($elementBedFile);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

print STDERR "running headers2bed ...\n" if($verbose);
my $headerBEDFile=headers2bed($matrixObject);
print STDERR "\t$headerBEDFile\n" if($verbose);

print STDERR "\n" if($verbose);

print STDERR "intersecting BED files ...\n" if($verbose);
my $elementFileName=getFileName($elementBedFile);
my $bedOverlapFile=intersectBED($headerBEDFile,$elementBedFile);
print STDERR "\t$bedOverlapFile\n" if($verbose);
system("rm '".$headerBEDFile."'");

print STDERR "\n" if($verbose);

print STDERR "loading BED file ...\n" if($verbose);
my ($elements)=loadBED($bedOverlapFile);
print STDERR "\tfound ".@{$elements}." elements\n" if($verbose);
system("rm '".$bedOverlapFile."'");

croak "found no overlapping headers!" if(@{$elements} == 0);

print STDERR "\n" if($verbose);

my $headerSizing=$matrixObject->{ headerSizing };

my $matrix={};
($matrix,$matrixObject)=getData($inputMatrix,$matrixObject,$verbose,$minDistance,(($maxDistance*2)+$headerSizing));

print STDERR "\n" if($verbose);

# calculate the boundaryReach index for each bin and store in a new data struct.
print STDERR "piling up data ...\n" if($verbose);
elementPileUp($matrixObject,$matrix,$elementFileName,$elements,$elementZoneSize,$minDistance,$maxDistance,$aggregrateMode,$excludeZero,$commentLine);

print STDERR "\n" if($verbose);
