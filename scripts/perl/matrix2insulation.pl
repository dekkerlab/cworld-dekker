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

    my ($inputMatrix,$verbose,$output,$insulationSquare,$smoothSize,$insulationDeltaSpan,$insulationMode,$noiseThreshold,$boundaryMarginOfError,$excludeZero,$yBound,$transparentBGFlag,$minDistance,$maxDistance,$meanInsulationScore);
    
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
    
    if( exists($opts->{ insulationSquare }) ) {
        $insulationSquare = $opts->{ insulationSquare };
        $insulationSquare=500000 if($insulationSquare <= 0);
    } else {
        $insulationSquare=500000;
    }
    
    if( exists($opts->{ smoothSize }) ) {
        $smoothSize = $opts->{ smoothSize };
    } else {
        $smoothSize=0;
    }
    
    if( exists($opts->{ insulationDeltaSpan }) ) {
        $insulationDeltaSpan = $opts->{ insulationDeltaSpan };
        $insulationDeltaSpan=($insulationSquare/2) if($insulationDeltaSpan <= 0);
    } else {
        $insulationDeltaSpan=($insulationSquare/2);
    }
    
    if( exists($opts->{ insulationMode }) ) {
        $insulationMode = $opts->{ insulationMode };
    } else {
        $insulationMode="mean";
    }
    
    if( exists($opts->{ noiseThreshold }) ) {
        $noiseThreshold = $opts->{ noiseThreshold };
    } else {
        $noiseThreshold=0;
    }
    
    if( exists($opts->{ boundaryMarginOfError }) ) {
        $boundaryMarginOfError = $opts->{ boundaryMarginOfError };
    } else {
        $boundaryMarginOfError=0;
    }
    
    if( exists($opts->{ excludeZero }) ) {
        $excludeZero = 1;
    } else {
        $excludeZero = 0;
    }
    
    if( exists($opts->{ yBound }) ) {
        $yBound = $opts->{ yBound };
    } else {
        $yBound = 0;
    }
    
    if( exists($opts->{ transparentBGFlag }) ) {
        $transparentBGFlag = 1;
    } else {
        $transparentBGFlag = 0;
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
    
    if( exists($opts->{ meanInsulationScore }) ) {
        $meanInsulationScore = $opts->{ meanInsulationScore };
    } else {
        $meanInsulationScore = "NA";
    }
    
    $ret->{ inputMatrix }=$inputMatrix;
    $ret->{ verbose }=$verbose;
    $ret->{ output }=$output;
    $ret->{ insulationSquare }=$insulationSquare;
    $ret->{ smoothSize }=$smoothSize;
    $ret->{ insulationDeltaSpan }=$insulationDeltaSpan;
    $ret->{ insulationMode }=$insulationMode;
    $ret->{ noiseThreshold }=$noiseThreshold;
    $ret->{ boundaryMarginOfError }=$boundaryMarginOfError;
    $ret->{ excludeZero }=$excludeZero;
    $ret->{ yBound }=$yBound;
    $ret->{ transparentBGFlag }=$transparentBGFlag;
    $ret->{ minDistance }=$minDistance;
    $ret->{ maxDistance }=$maxDistance;
    $ret->{ meanInsulationScore }=$meanInsulationScore;
    
    return($ret,$inputMatrix,$verbose,$output,$insulationSquare,$smoothSize,$insulationDeltaSpan,$insulationMode,$noiseThreshold,$boundaryMarginOfError,$excludeZero,$yBound,$transparentBGFlag,$minDistance,$maxDistance,$meanInsulationScore);
}

sub intro() {
    print STDERR "\n";
    
    print STDERR "Tool:\t\t".$tool."\n";
    print STDERR "Version:\t".$cworld::dekker::VERSION."\n";
    print STDERR "Summary:\tcalculate insulation index (TADs) of supplied matrix\n";
    
    print STDERR "\n";
}

sub help() {
    intro();
    
    print STDERR "Usage: perl matrix2insulation.pl [OPTIONS] -i <inputMatrix>\n";
    
    print STDERR "\n";
    
    print STDERR "Required:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-i", "[]", "input matrix file");
    
    print STDERR "\n";
    
    print STDERR "Options:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-v", "[]", "FLAG, verbose mode");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-o", "[]", "prefix for output file(s)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--is", "[500000]", "optional, insulation square size, size of the insulation square in bp");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ss", "[0]", "optional, insulation smooth size, size of insulation vector smoothing vector in bp");   
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ids", "[250000]", "optional, insulation delta span (window size of insulationd delta - blue line)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--im", "[mean]", "optional, insulation mode (how to aggregrate signal within insulation square), mean,sum,median");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--nt", "[0.1]", "optional, noise threshold, minimum depth of valley");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--bmoe", "[3]", "optional, boundary margin of error (specified in number of BINS), added to each side of the boundary");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--yb", "[auto]", "optional, -yBound - +yBound of insulation plot");    
    printf STDERR ("\t%-10s %-10s %-10s\n", "--bg", "[]", "FLAG, use transparent background of insulation plot");    
    printf STDERR ("\t%-10s %-10s %-10s\n", "--minDist", "[]", "minimum allowed interaction distance, exclude < N distance (in BP)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--maxDist", "[]", "maximum allowed interaction distance, exclude > N distance (in BP)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ez", "[]", "FLAG, ignore 0s in all calculations");    
   
    print STDERR "\n";
    
    print STDERR "Notes:";
    print STDERR "
    This script calculates the insulation index of a given matrix to identify TAD boundaries.
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
    my $insulationFile=shift;
    my $insulationSquare=shift;
    my $smoothSize_binSize=shift;
    my $insulationMode=shift;
    my $excludeZero=shift;
    my $meanInsulationScore=shift;
    
    my $inc2header=$matrixObject->{ inc2header };
    my $header2inc=$matrixObject->{ header2inc };
    my $numHeaders=$matrixObject->{ numYHeaders };
    my $headerSizing=$matrixObject->{ headerSizing };
    my $headerSpacing=$matrixObject->{ headerSpacing };
    my $verbose=$matrixObject->{ verbose };
    
    my %matrixInsulation=();
    my @insulationSignals=();
    
    my $insulationStartIndex=$insulationSquare;
    my $insulationEndIndex=($numHeaders-$insulationSquare);
    
    print STDERR "\traw insulation\n" if($verbose);
    for(my $y=$insulationStartIndex;$y<$insulationEndIndex;$y++) {
        my $yHead=$inc2header->{ y }->{$y};
        
        #do not include the diagonal
        my $startY=($y-$insulationSquare);
        $startY=0 if($startY < 0);
        my $endY=$y;        
        my $startX=$y+1;
        $startX=$numHeaders if($startX > $numHeaders);
        my $endX=($y+$insulationSquare)+1;
        $endX=$numHeaders if($endX > $numHeaders);
            
        my $skipFlag=0;
        my @boxData=();
        
        my $expectedDataCount=(($endY-$startY)*($endX-$startX));
        my $nullCount=0;
        my $zeroCount=0;
        
        for(my $y2=$startY;$y2<$endY;$y2++) {
            my $yHead=$inc2header->{ y }->{$y2};
            for(my $x2=$startX;$x2<$endX;$x2++) {    
                my $xHead=$inc2header->{ x }->{$x2};
                
                my $inten=$matrixObject->{ missingValue };
                $inten=$matrix->{$y2}->{$x2} if(exists($matrix->{$y2}->{$x2}));
                
                $nullCount++ if($inten eq "NA");
                $zeroCount++ if($inten eq 0);
                
                $skipFlag = 1 if($inten eq "NA");
                next if($inten eq "NA");
                
                next if(($inten eq 0) and ($excludeZero));
                
                push(@boxData,$inten);
            }
        }
        
        my $boxDataSize=scalar @boxData;
        
        my $maxMissing=($expectedDataCount/2);
        my $maxZero=($expectedDataCount/2);
        
        next if($nullCount > $maxMissing); # skip if high NAs
        next if(($excludeZero) and ($zeroCount > $maxMissing)); # skip if high 0s
        
        next if(($skipFlag) and ($insulationMode eq "sum"));
        next if($boxDataSize == 0);
        
        my $boxDataStats=listStats(\@boxData);
        my $boxDataAggregate="NA";
        $boxDataAggregate=$boxDataStats->{ $insulationMode } if(exists($boxDataStats->{ $insulationMode }));
        
        # store result in hash header->insulation
        $matrixInsulation{$yHead}{ insulation }=$boxDataAggregate;
        push(@insulationSignals,$boxDataAggregate) if($boxDataAggregate ne "NA");
    }
    
    my $smoothWindow=$smoothSize_binSize;
    print STDERR "\tsmoothed insulation ($smoothWindow)\n" if($verbose);
    for(my $y=$insulationStartIndex;$y<$insulationEndIndex;$y++) {
        my $yHead=$inc2header->{ y }->{$y};
        
        next if(!exists($matrixInsulation{$yHead}{ insulation }));
        
        my $insulationScore=$matrixInsulation{$yHead}{ insulation };
        
        my @tmp=();
        for(my $i=($y-$smoothWindow);$i<=($y+$smoothWindow);$i++) {
            next if( ($i < 0) or ($i > $insulationEndIndex) );
            
            my $tmpHead=$inc2header->{ y }->{$i};
            next if(!exists($matrixInsulation{$tmpHead}{ insulation }));
            my $tmpScore=$matrixInsulation{$tmpHead}{ insulation };
            push(@tmp,$tmpScore) if($tmpScore ne "NA");
        }
        my $tmpStats=listStats(\@tmp) if(@tmp > 0);
        my $smoothedScore="NA";
        $smoothedScore=$tmpStats->{ mean } if(exists($tmpStats->{ mean }));
        $matrixInsulation{$yHead}{ smoothed_insulation }=$smoothedScore;
    }
    
    print STDERR "\tnormalized insulation\n" if($verbose);
    
    # mean center the insulation data
    croak "no available data points" if(@insulationSignals == 0);
    
    my $tmpStats=listStats(\@insulationSignals);
    $meanInsulationScore=$tmpStats->{ mean } if($meanInsulationScore eq "NA");
    print STDERR "\tmean centering data ($meanInsulationScore)\n" if($verbose);
    
    my @normalizedInsulationSignals=();
    for(my $y=$insulationStartIndex;$y<$insulationEndIndex;$y++) {
        my $yHead=$inc2header->{ y }->{$y};
        
        next if(!exists($matrixInsulation{$yHead}{ smoothed_insulation }));
        
        my $insulationScore=$matrixInsulation{$yHead}{ smoothed_insulation };
        
        #mean center the insulation scores
        my $normalizedInsulationScore = "NA";
        $normalizedInsulationScore = (log($insulationScore/$meanInsulationScore)/log(2)) if(($meanInsulationScore ne "NA") and ($insulationScore ne "NA") and ($meanInsulationScore != 0) and ($insulationScore != 0));
         
        $matrixInsulation{$yHead}{ normalized_insulation }=$normalizedInsulationScore;
        push(@normalizedInsulationSignals,$normalizedInsulationScore) if($normalizedInsulationScore ne "NA");
    }
    
    my $tmpNormalizedStats=listStats(\@normalizedInsulationSignals);
    my $meanNormalizedInsulationScore=$tmpNormalizedStats->{ mean };
    print STDERR "\tnormalized data mean ($meanNormalizedInsulationScore)\n" if($verbose);
    
    return(\%matrixInsulation);
    
}

sub outputInsulation($$$$$$$) {
    my $matrixObject=shift;
    my $matrixInsulation=shift;
    my $matrixDelta=shift;
    my $matrixDeltaSquare=shift;    
    my $insulationFile=shift;
    my $yBound=shift;
    my $commentLine=shift;
    
    my $inc2header=$matrixObject->{ inc2header };
    my $header2inc=$matrixObject->{ header2inc };
    my $numHeaders=$matrixObject->{ numYHeaders };
    my $headerSizing=$matrixObject->{ headerSizing };
    my $headerSpacing=$matrixObject->{ headerSpacing };
    
    # main insulation file
    open(OUT,outputWrapper($insulationFile,$commentLine)) or croak "Could not open file [$insulationFile] - $!";
    print OUT "header\tstart\tend\tmidpoint\tbinStart\tbinEnd\tbinMidpoint\trawInsulationScore\tsmoothedInsulaton\tinsulationScore\tdelta\tdeltaSquare\n";
    
    # bed graph file of insulation data
    open(BEDGRAPH,outputWrapper($insulationFile.".bedGraph")) or croak "Could not open file [$insulationFile] - $!";
    print BEDGRAPH "track type=bedGraph name='".$insulationFile."' description='".$insulationFile." - insutation score' maxHeightPixels=128:64:32 visibility=full autoScale=off viewLimits=-".$yBound.":".$yBound." color=0,0,0 altColor=100,100,100\n" if($yBound != 0);
    print BEDGRAPH "track type=bedGraph name='".$insulationFile."' description='".$insulationFile." - insutation score' maxHeightPixels=128:64:32 visibility=full autoScale=on color=0,0,0 altColor=100,100,100\n" if($yBound == 0);
    
    for(my $y=0;$y<$numHeaders;$y++) {
        # dump insulation data to file
        my $yHead=$inc2header->{ y }->{$y};
        
        my $insulation="NA";
        $insulation=$matrixInsulation->{$yHead}{ insulation } if(exists($matrixInsulation->{$yHead}{ insulation }));
        
        my $smoothedInsulation="NA";
        $smoothedInsulation=$matrixInsulation->{$yHead}{ smoothed_insulation } if(exists($matrixInsulation->{$yHead}{ smoothed_insulation }));
        
        my $normalizedInsulation="NA";
        $normalizedInsulation=$matrixInsulation->{$yHead}{ normalized_insulation } if(exists($matrixInsulation->{$yHead}{ normalized_insulation }));
        
        my $delta="NA";
        $delta=$matrixDelta->{$yHead} if(exists($matrixDelta->{$yHead}));
        
        my $deltaSquare="NA";
        $deltaSquare=$matrixDeltaSquare->{$yHead} if(exists($matrixDeltaSquare->{$yHead}));
        
        my $yHeadObject=getHeaderObject($yHead);
        my $yHeadChromosome=$yHeadObject->{ chromosome };
        my $yHeadStart=$yHeadObject->{ start };
        my $yHeadEnd=$yHeadObject->{ end };
        my $yHeadMidpoint=round(($yHeadStart+$yHeadEnd)/2);
        
        my $binStart = round($yHeadStart/$headerSpacing);
        my $binEnd = round($yHeadEnd/$headerSpacing);
        my $binMidpoint=(($binStart+$binEnd)/2);
        
        print OUT "$yHead\t$yHeadStart\t$yHeadEnd\t$yHeadMidpoint\t$binStart\t$binEnd\t$binMidpoint\t$insulation\t$smoothedInsulation\t$normalizedInsulation\t$delta\t$deltaSquare\n";
        
        # strip off chr group if exists for proper UCSC usage
        $yHeadChromosome=stripChromosomeGroup($yHeadChromosome);
        
        $insulation=0 if($insulation eq "NA");
        print BEDGRAPH "$yHeadChromosome\t$yHeadStart\t$yHeadEnd\t$normalizedInsulation\n";
    }
    
    close(OUT);
    close(BEDGRAPH);
    
}

sub outputInsulationBoundaries($$$$) {
    my $matrixObject=shift;
    my $tadBoundaries=shift;
    my $insulationFile=shift;
    my $commentLine=shift;
    
    my $inc2header=$matrixObject->{ inc2header };
    my $header2inc=$matrixObject->{ header2inc };
    my $numHeaders=$matrixObject->{ numYHeaders };
    my $headerSizing=$matrixObject->{ headerSizing };
    my $headerSpacing=$matrixObject->{ headerSpacing };

    open(OUT,outputWrapper($insulationFile,$commentLine)) or croak "Could not open file [$insulationFile] - $!";
    print OUT "header\tstart\tend\tbinStart\tbinEnd\tbinMidpoint\tboundaryHeader\tboundaryStrength\tboundaryDirectionality\tboundaryInsulation\n";
    
    open(BED,outputWrapper($insulationFile.".bed")) or croak "Could not open file [$insulationFile] - $!";
    
    print BED "track name='".$insulationFile."' description='called tad boundaries'\n";
    
    for(my $y=0;$y<$numHeaders;$y++) {
        # dump insulation data to file
        my $yHead=$inc2header->{ y }->{$y};
    
        next if(!exists($tadBoundaries->{$yHead}->{ header }));
        
        my $boundaryHeader=$tadBoundaries->{$yHead}->{ header };
        next if($boundaryHeader eq "NA");
        
        my $boundaryStrength=$tadBoundaries->{$yHead}->{ strength };
        my $boundaryInsulation=$tadBoundaries->{$yHead}->{ insulation };
        my $boundaryDirectionality=$tadBoundaries->{$yHead}->{ directionality };
        
        my $boundaryObject=getHeaderObject($boundaryHeader);
        my $boundaryChromosome=$boundaryObject->{ chromosome };
        my $boundaryStart=$boundaryObject->{ start };
        my $boundaryEnd=$boundaryObject->{ end };
        $boundaryChromosome=stripChromosomeGroup($boundaryChromosome);
        
        print BED "$boundaryChromosome\t$boundaryStart\t$boundaryEnd\t$yHead\t$boundaryStrength\n";

        my $binStart = round($boundaryStart/$headerSpacing);
        my $binEnd = round($boundaryEnd/$headerSpacing);
        my $binMidpoint=(($binStart+$binEnd)/2);
        
        print OUT "$boundaryHeader\t$boundaryStart\t$boundaryEnd\t$binStart\t$binEnd\t$binMidpoint\t$yHead\t$boundaryStrength\t$boundaryDirectionality\t$boundaryInsulation\n";
        
    }
    
    close(OUT);
    close(BED);
}

sub detectInsulationBoundaries($$$$$$$) {
    my $matrixObject=shift;
    my $inputMatrixName=shift;
    my $matrixInsulation=shift;
    my $insulationSquare=shift;
    my $insulationDeltaSpan_binSize=shift;
    my $noiseThreshold=shift;
    my $boundaryMarginOfError=shift;
    
    my $inc2header=$matrixObject->{ inc2header };
    my $header2inc=$matrixObject->{ header2inc };
    my $numHeaders=$matrixObject->{ numYHeaders };
    my $headerSizing=$matrixObject->{ headerSizing };
    my $headerSpacing=$matrixObject->{ headerSpacing };
    my $output=$matrixObject->{ output };
    
    my $halfInsulationDeltaSpan_binSize=ceil($insulationDeltaSpan_binSize/2);
    
    my %tadBoundaries=();
    my $matrixDelta={};
    
    for(my $y=0;$y<$numHeaders;$y++) {
        # dump insulation data to file
        my $yHead=$inc2header->{ y }->{$y};
        
        my $normalized_insulation="NA";
        $normalized_insulation=$matrixInsulation->{$yHead}{ normalized_insulation } if(exists($matrixInsulation->{$yHead}{ normalized_insulation }));
        
        my $leftBound=($y-$halfInsulationDeltaSpan_binSize);
        $leftBound=0 if($leftBound < 0);
        my $rightBound=($y+$halfInsulationDeltaSpan_binSize);
        $rightBound=$numHeaders-1 if($rightBound >= $numHeaders);
                
        my @leftDelta=();
        for(my $ly=$leftBound;$ly<$y;$ly++) {
            my $leftHead=$inc2header->{ y }->{$ly};
            my $leftInsulation="NA";
            $leftInsulation=$matrixInsulation->{$leftHead}{ normalized_insulation } if(exists($matrixInsulation->{$leftHead}{ normalized_insulation }));
            
            my $delta="NA";
            $delta=($leftInsulation-$normalized_insulation) if(($normalized_insulation ne "NA") and ($leftInsulation ne "NA"));
            push(@leftDelta,$delta) if($delta ne "NA");
        }
        my $leftDeltaAggregrate="NA";
        if(@leftDelta > 0) {
            my $leftDeltaStats=listStats(\@leftDelta);
            $leftDeltaAggregrate=$leftDeltaStats->{ mean };
        }
        
        my @rightDelta=();
        for(my $ry=$y+1;$ry<=$rightBound;$ry++) {
            my $rightHead=$inc2header->{ y }->{$ry};
            my $rightInsulation="NA";
            $rightInsulation=$matrixInsulation->{$rightHead}{ normalized_insulation } if(exists($matrixInsulation->{$rightHead}{ normalized_insulation }));
            
            my $delta="NA";
            $delta=($rightInsulation-$normalized_insulation) if(($normalized_insulation ne "NA") and ($rightInsulation ne "NA"));
            push(@rightDelta,$delta) if($delta ne "NA");
        }
        my $rightDeltaAggregrate="NA";
        if(@rightDelta > 0) {
            my $rightDeltaStats=listStats(\@rightDelta);
            $rightDeltaAggregrate=$rightDeltaStats->{ mean };
        }
        
        my $deltaDelta="NA";
        $deltaDelta=($leftDeltaAggregrate-$rightDeltaAggregrate) if(($leftDeltaAggregrate ne "NA") and ($rightDeltaAggregrate ne "NA"));        
        $matrixDelta->{$yHead}=$deltaDelta;
    }
    
    my %matrixDeltaSquare=();
    for(my $y=0;$y<$numHeaders;$y++) {
        # dump insulation data to file
        my $yHead=$inc2header->{ y }->{$y};
        my $delta="NA";
        $delta=$matrixDelta->{$yHead} if(exists($matrixDelta->{$yHead}));
        
        $matrixDeltaSquare{$yHead}=$delta;
        $matrixDeltaSquare{$yHead}=1 if(($delta ne "NA") and ($delta > 0));
        $matrixDeltaSquare{$yHead}=-1 if(($delta ne "NA") and ($delta < 0));
    }
    
    # filter out the peaks
    my @preLimBoundaries=();
    my $binc=0;
    my $lastBound=0;
    for(my $y=0;$y<$numHeaders;$y++) {
        # dump insulation data to file
        my $yHead=$inc2header->{ y }->{$y};
        my $deltaSquare="NA";
        $deltaSquare=$matrixDeltaSquare{$yHead} if(exists($matrixDeltaSquare{$yHead}));
        
        next if($deltaSquare eq "NA");
                
        if(($deltaSquare != $lastBound) and ($lastBound != 0)) {
            $preLimBoundaries[$binc]{ header }=$yHead;
            $preLimBoundaries[$binc]{ lastBound }=$lastBound;
            $binc++;
        } 
        $lastBound=$deltaSquare;
    }
    
    # for all valleys
    my @boundaryStrengthArr=();
    my @boundaryDirectionalityArr=();
    my $boundary_num=0;
    for(my $i=0;$i<@preLimBoundaries;$i++) {
        my $yHead=$preLimBoundaries[$i]{ header };
        my $yIndex=$header2inc->{ y }->{$yHead};
        
        my $currentDelta="NA";
        $currentDelta=$matrixDelta->{$yHead} if(exists($matrixDelta->{$yHead}));
        
        next if($currentDelta eq "NA");
        
        my ($leftSearchIndex,$rightSearchIndex);
        
        # get strength via delta
        
        $leftSearchIndex=$yIndex;
        # while left is increasing
        while( ($leftSearchIndex > 0) and ($matrixDelta->{$inc2header->{ y }->{$leftSearchIndex}} ne "NA") and ($matrixDelta->{$inc2header->{ y }->{$leftSearchIndex-1}} ne "NA") and ($matrixDelta->{$inc2header->{ y }->{$leftSearchIndex-1}} >= $matrixDelta->{$inc2header->{ y }->{$leftSearchIndex}}) ) {
            $leftSearchIndex--;
        }
        my $leftDeltaBound="NA";
        $leftDeltaBound=$matrixDelta->{$inc2header->{ y }->{$leftSearchIndex}} if(exists($matrixDelta->{$inc2header->{ y }->{$leftSearchIndex}}));
        
        $rightSearchIndex=$yIndex;
        # while right is decreasing
        while( ($rightSearchIndex < ($numHeaders-2)) and ($matrixDelta->{$inc2header->{ y }->{$rightSearchIndex}} ne "NA") and ($matrixDelta->{$inc2header->{ y }->{$rightSearchIndex+1}} ne "NA") and ($matrixDelta->{$inc2header->{ y }->{$rightSearchIndex+1}} <= $matrixDelta->{$inc2header->{ y }->{$rightSearchIndex}}) ) {
            $rightSearchIndex++;
        }
        my $rightDeltaBound="NA";
        $rightDeltaBound=$matrixDelta->{$inc2header->{ y }->{$rightSearchIndex}};
        
        my $valleyDeltaStrength="NA";
        $valleyDeltaStrength=($leftDeltaBound-$rightDeltaBound) if(($leftDeltaBound ne "NA") and ($rightDeltaBound ne "NA"));
                
        next if($valleyDeltaStrength <= 0); # remove peaks
        next if($valleyDeltaStrength eq "NA");
        next if($valleyDeltaStrength < $noiseThreshold);
                
        my $valleyDirectionality="NA";
        $valleyDirectionality=log((abs($leftDeltaBound)/abs($rightDeltaBound)))/log(2) if(($leftDeltaBound != 0) and ($rightDeltaBound != 0));
        
        my $boundaryObject=getHeaderObject($yHead);
        my $boundaryChromosome=$boundaryObject->{ chromosome };
        my $boundaryAssembly=$boundaryObject->{ assembly };
        my $boundaryStart=$boundaryObject->{ start };
        my $boundaryEnd=$boundaryObject->{ end };
        my $boundaryMidpoint=round(($boundaryStart+$boundaryEnd)/2);
        
        $boundaryStart -= ($boundaryMarginOfError*$headerSpacing);
        $boundaryEnd += ($boundaryMarginOfError*$headerSpacing);
        
        my $normalized_insulation="NA";
        $normalized_insulation=$matrixInsulation->{$yHead}{ normalized_insulation } if(exists($matrixInsulation->{$yHead}{ normalized_insulation }));
        
        my $boundaryHeader="boundary.".$boundary_num."|".$boundaryAssembly."|".$boundaryChromosome.":".$boundaryStart."-".$boundaryEnd;
        $tadBoundaries{$yHead}{ header }=$boundaryHeader;
        $tadBoundaries{$yHead}{ strength }=$valleyDeltaStrength;
        $tadBoundaries{$yHead}{ insulation }=$normalized_insulation;
        $tadBoundaries{$yHead}{ directionality }=$valleyDirectionality;
        
        push(@boundaryStrengthArr,$valleyDeltaStrength) if($valleyDeltaStrength ne "NA");
        push(@boundaryDirectionalityArr,$valleyDirectionality) if($valleyDirectionality ne "NA");
        $boundary_num++;
    }
    
    # log file
    my $logFile=$output.".log";
    open(LOG,outputWrapper($logFile)) or croak "Could not open file [$logFile] - $!";
    
    print LOG "insulationSquare\tnumBoundaries\tboundaryStrength_mean\tboundaryStrength_median\tboundaryStrength_spread\tboundaryStrength_stdev\tboundaryStrength_iqr\tboundaryStrength_iqrMean\tboundaryStrength_iqrStdev\tboundaryStrength_iqrVariance\tboundaryDirectionality_mean\n";
    
    my $numBoundaries=keys %tadBoundaries;
    
    my $boundaryStrengthStats=listStats(\@boundaryStrengthArr);
    my $boundaryStrength_mean=$boundaryStrengthStats->{ mean };
    my $boundaryStrength_median=$boundaryStrengthStats->{ median };
    my $boundaryStrength_spread=$boundaryStrengthStats->{ variance };
    my $boundaryStrength_stdev=$boundaryStrengthStats->{ stdev };
    my $boundaryStrength_iqr=$boundaryStrengthStats->{ iqr };
    my $boundaryStrength_iqrStdev=$boundaryStrengthStats->{ iqrStdev };
    my $boundaryStrength_iqrVariance=$boundaryStrengthStats->{ iqrVariance };
    my $boundaryStrength_iqrMean=$boundaryStrengthStats->{ iqrMean };
    
    my $boundaryDirectionalityStats=listStats(\@boundaryDirectionalityArr);
    my $boundaryDirectionality_mean=$boundaryDirectionalityStats->{ mean };
    
    print LOG "$insulationSquare\t$numBoundaries\t$boundaryStrength_mean\t$boundaryStrength_median\t$boundaryStrength_spread\t$boundaryStrength_stdev\t$boundaryStrength_iqr\t$boundaryStrength_iqrMean\t$boundaryStrength_iqrStdev\t$boundaryStrength_iqrVariance\t$boundaryDirectionality_mean\n";
    
    close(LOG);
    
    return(\%tadBoundaries,$matrixDelta,\%matrixDeltaSquare);
}

my %options;
my $results = GetOptions( \%options,'inputMatrix|i=s','verbose|v','output|o=s','insulationSquare|is=i','smoothSize|ss=i','insulationDeltaSpan|ids=i','insulationMode|im=s','noiseThreshold|nt=f','boundaryMarginOfError|bmoe=i','excludeZero|ez','yBound|yb=f','transparentBGFlag|bg','minDistance|minDist=i','maxDistance|maxDist=i','meanInsulationScore|mis=f') or croak help();
my ($ret,$inputMatrix,$verbose,$output,$insulationSquare,$smoothSize,$insulationDeltaSpan,$insulationMode,$noiseThreshold,$boundaryMarginOfError,$excludeZero,$yBound,$transparentBGFlag,$minDistance,$maxDistance,$meanInsulationScore)=check_options( \%options );

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

croak "matrix must be symmetrical" if(!$symmetrical);
croak "matrix must be equally-binned" if(!$equalHeaderFlag);

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
    
print STDERR "\n" if($verbose);

# get fragment spacing (i.e. bin size)
my ($equalSpacingFlag,$equalSizingFlag,$headerSpacing,$headerSizing)=getHeaderSpacing($inc2header->{ y });

# insulation square size
my $insulationSquare_binSize=ceil(($insulationSquare-($headerSizing-$headerSpacing))/$headerSpacing);
$insulationSquare_binSize = 2 if($insulationSquare_binSize <= 1);
my $insulationSquare_bpSize=($insulationSquare_binSize * $headerSpacing)+($headerSizing-$headerSpacing);
print STDERR "insulationSquare\t$insulationSquare\t$insulationSquare_binSize\t$insulationSquare_bpSize\n" if($verbose);
croak "insulationSquare_binSize cannot be larger than dataset ($insulationSquare_binSize > $numHeaders)" if($insulationSquare_binSize > $numHeaders);

# set delta size
my $insulationDeltaSpan_binSize=ceil(($insulationDeltaSpan)/$headerSpacing);
$insulationDeltaSpan_binSize += 1 if(($insulationDeltaSpan_binSize % 2));
my $insulationDeltaSpan_bpSize=($insulationDeltaSpan_binSize * $headerSpacing)+($headerSizing-$headerSpacing);
print STDERR "insulationDeltaSpan\t$insulationDeltaSpan\t$insulationDeltaSpan_binSize\t$insulationDeltaSpan_bpSize\n" if($verbose);
croak "insulationDeltaSpan cannot be larger than dataset ($insulationDeltaSpan_binSize > $numHeaders)" if($insulationDeltaSpan_binSize > $numHeaders);

# set smooth size
my $smoothSize_binSize=ceil(($smoothSize-($headerSizing-$headerSpacing))/$headerSpacing);
$smoothSize_binSize=0 if($smoothSize_binSize < 0);
my $smoothSize_bpSize=($smoothSize_binSize * $headerSpacing)+($headerSizing-$headerSpacing);
print STDERR "smoothSize\t$smoothSize\t$smoothSize_binSize\t$smoothSize_bpSize\n" if($verbose);
croak "smoothSize_binSize cannot be larger than dataset ($smoothSize_binSize > $numHeaders)" if($smoothSize_binSize > $numHeaders);

print STDERR "\n" if($verbose);

# add suffix
$output .= "--is".$insulationSquare_bpSize."--nt".$noiseThreshold."--ids".$insulationDeltaSpan_bpSize."--ss".$smoothSize_bpSize."--im".$insulationMode;

carp "insulation image size is < numHeaders ($numHeaders), visual artifacts are possible!" if($numHeaders > $imageWidth);

#read Matrix
my ($matrix)=getData($inputMatrix,$matrixObject,$verbose,$minDistance,(($insulationSquare_bpSize*2)+$headerSizing));

# reset sparse value to NA
$matrixObject->{ missingValue }="NA";
$missingValue=$matrixObject->{ missingValue };

print STDERR "\n" if($verbose);

# calculate the insulation index for each bin and store in a new data struct.
print STDERR "calculating insulation index...\n" if($verbose);
my ($matrixInsulation)=calculateInsulation($matrixObject,$matrix,$inputMatrixName,$insulationSquare_binSize,$smoothSize_binSize,$insulationMode,$excludeZero,$meanInsulationScore);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

# detect boundaries
print STDERR "detecting boundaries...\n" if($verbose);
my ($tadBoundaries,$matrixDelta,$matrixDeltaSquare)=detectInsulationBoundaries($matrixObject,$inputMatrixName,$matrixInsulation,$insulationSquare_bpSize,$insulationDeltaSpan_binSize,$noiseThreshold,$boundaryMarginOfError);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

# write insulation data to file
my $insulationFile=$output.".insulation";
print STDERR "outputing insulation...\n" if($verbose);
outputInsulation($matrixObject,$matrixInsulation,$matrixDelta,$matrixDeltaSquare,$insulationFile,$yBound,$commentLine);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

# write insulation data to file
print STDERR "outputing insulation boundaries...\n" if($verbose);
outputInsulationBoundaries($matrixObject,$tadBoundaries,$insulationFile.".boundaries",$commentLine);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

if($numContigs == 1) {
    # plot insulation using R
    print STDERR "plotting insulation...\n" if($verbose);
    system("Rscript '".$scriptPath."/R/matrix2insulation.R' '".$cwd."' '".$insulationFile."' ".$headerSizing." ".$headerSpacing." ".$insulationSquare." ".$insulationSquare_bpSize." ".$imageWidth." ".$yBound." ".$transparentBGFlag." > /dev/null");
    print STDERR "\tdone\n" if($verbose);
}

print STDERR "\n" if($verbose);