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

    my ($inputMatrix,$verbose,$output,$elementBedFiles,$interactionBedFile,$elementName,$elementZoneSize,$minDistance,$maxDistance,$minElementDistance,$maxElementDistance,$aggregrateMode,$excludeZero,$includeDiagonal,$excludeCis,$excludeTrans,$debugMode);
    
    my $ret={};
    
    if( exists($opts->{ inputMatrix }) ) {
        $inputMatrix = $opts->{ inputMatrix };
    } else {
        print STDERR "\nERROR: Option inputMatrix|i is required.\n";
        help();
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
    
    if( exists($opts->{ elementBedFiles }) ) {
        $elementBedFiles = $opts->{ elementBedFiles };
    } else {
        $elementBedFiles=[];
    }
    
       if( exists($opts->{ interactionBedFile }) ) {
        $interactionBedFile = $opts->{ interactionBedFile };
    } else {
        $interactionBedFile="";
    }
    
    if( exists($opts->{ elementName }) ) {
        $elementName = $opts->{ elementName };
    } else {
        $elementName=getSmallUniqueString();
    }
    
    if( exists($opts->{ elementZoneSize }) ) {
        $elementZoneSize = $opts->{ elementZoneSize };
    } else {
        $elementZoneSize=100000;
    }
    
    if( exists($opts->{ minDistance }) ) {
        $minDistance = $opts->{ minDistance };
    } else {
        $minDistance=undef;
    }
    
    if( exists($opts->{ maxDistance }) ) {
        $maxDistance = $opts->{ maxDistance };
    } else {
        $maxDistance=undef;
    }
    
    if( exists($opts->{ minElementDistance }) ) {
        $minElementDistance = $opts->{ minElementDistance };
    } else {
        $minElementDistance=undef;
    }
    
    if( exists($opts->{ maxElementDistance }) ) {
        $maxElementDistance = $opts->{ maxElementDistance };
    } else {
        $maxElementDistance=undef;
    }
    
    if( exists($opts->{ aggregrateMode }) ) {
        $aggregrateMode = $opts->{ aggregrateMode };
    } else {
        $aggregrateMode="mean";
    }
    
    if( exists($opts->{ excludeZero }) ) {
        $excludeZero = 1;
    } else {
        $excludeZero = 0;
    }
    
    if( exists($opts->{ includeDiagonal }) ) {
        $includeDiagonal = 1;
    } else {
        $includeDiagonal = 0;
    }
    
    if( exists($opts->{ excludeCis }) ) {
        $excludeCis = 1;
    } else {
        $excludeCis = 0;
    }
    
    if( exists($opts->{ excludeTrans }) ) {
        $excludeTrans = 1;
    } else {
        $excludeTrans = 0;
    }
    
    if( exists($opts->{ debugMode }) ) {
        $debugMode = 1;
    } else {
        $debugMode = 0;
    }
    
    $ret->{ inputMatrix }=$inputMatrix;
    $ret->{ verbose }=$verbose;
    $ret->{ output }=$output;
    $ret->{ elementBedFiles }=$elementBedFiles;
    $ret->{ interactionBedFile }=$interactionBedFile;
    $ret->{ elementName }=$elementName;
    $ret->{ elementZoneSize }=$elementZoneSize;
    $ret->{ minDistance }=$minDistance;
    $ret->{ maxDistance }=$maxDistance;
    $ret->{ minElementDistance }=$minElementDistance;
    $ret->{ maxElementDistance }=$maxElementDistance;
    $ret->{ aggregrateMode }=$aggregrateMode;
    $ret->{ excludeZero }=$excludeZero;
    $ret->{ includeDiagonal }=$includeDiagonal;
    $ret->{ excludeCis }=$excludeCis;
    $ret->{ excludeTrans }=$excludeTrans;
    $ret->{ debugMode }=$debugMode;
    
    return($ret,$inputMatrix,$verbose,$output,$elementBedFiles,$interactionBedFile,$elementName,$elementZoneSize,$minDistance,$maxDistance,$minElementDistance,$maxElementDistance,$aggregrateMode,$excludeZero,$includeDiagonal,$excludeCis,$excludeTrans,$debugMode);
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
    
    print STDERR "Usage: perl interactionPileUp.pl [OPTIONS] -i <inputMatrix> -ebf <elementBedFile>\n";
    
    print STDERR "\n";
    
    print STDERR "Required:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-i", "[]", "input matrix file");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ebf@", "[]", "element bed file (can accept multiple files)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ibf@", "[]", "interaction file, 7 column, 2 x bed3 + name");
    
    print STDERR "\n";
    
    print STDERR "Options:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-v", "[]", "FLAG, verbose mode");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-o", "[]", "prefix for output file(s)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--en", "[]", "elementName, descriptor for output files");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ezs", "[100000]", "elementZoneSize, size of zone to use around element (x axis - in BP)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--am", "[]", "aggregrateMode, how to aggregrate data (mean,sum,median,iqrMean,min,max etc)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--minDist", "[]", "minimum allowed interaction distance, exclude < N distance (in BP)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--maxDist", "[]", "maximum allowed interaction distance, exclude > N distance (in BP)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--minED", "[]", "minElementDistance, minimum element distance (#elements) to consider");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--maxED", "[]", "maxElementDistance, maximum element distance (#elements) to consider");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ez", "[]", "FLAG, exclude 0s in all calculations");    
    printf STDERR ("\t%-10s %-10s %-10s\n", "--id", "[off]", "FLAG, include any interaction space that touches diagonal");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ec", "[]", "FLAG, exclude CIS interactions");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--et", "[]", "FLAG, exclude TRANS interactions");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--d", "[]", "FLAG, debug mode (print out all element-element matrices)");

    print STDERR "\n";
    
    print STDERR "Notes:";
    print STDERR "
    This script aggregrates cData surrounding a list of 'elements'.
    Input Matrix can be TXT or gzipped TXT.
    AnchorListFile must be Bed5+ format.
    See website for matrix format details.\n";
    
    print STDERR "\n";
    
    print STDERR "Contact:
    Dekker Lab
    http://my5C.umassmed.edu
    my5C.help\@umassmed.edu\n";
    
    print STDERR "\n";
    
    exit;
}

sub interactionPileUp($$$$$$$$$$$$$$$$$$) {
    my $matrixObject=shift;
    my $matrix=shift;
    my $elementFileName=shift;
    my $elements=shift;
    my $interactions=shift;
    my $elementZoneSize=shift;
    my $minDistance=shift;
    my $maxDistance=shift;
    my $minElementDistance=shift;
    my $maxElementDistance=shift;
    my $aggregrateMode=shift;
    my $excludeZero=shift;
    my $includeDiagonal=shift;
    my $excludeCis=shift;
    my $excludeTrans=shift;
    my $commentLine=shift;
    my $debugMode=shift;
    my $verbose=shift;
    
    my $inc2header=$matrixObject->{ inc2header };
    my $header2inc=$matrixObject->{ header2inc };
    my $numXHeaders=$matrixObject->{ numXHeaders };
    my $numYHeaders=$matrixObject->{ numYHeaders };
    my $numTotalHeaders=$matrixObject->{ numTotalHeaders };
    my $headerSizing=$matrixObject->{ headerSizing };
    my $headerSpacing=$matrixObject->{ headerSpacing };
    my $output=$matrixObject->{ output };
    
    my $maxDistanceLimit=($numTotalHeaders * $headerSpacing)+($headerSizing-$headerSpacing);
    my $maxDistanceLimit_bins=$numTotalHeaders;
    
    # do not allow either element/distance to reach beyond matrix bounds
    $elementZoneSize=$maxDistanceLimit if($elementZoneSize > $maxDistanceLimit);
    
    my $elementZoneSize_bins=ceil(($elementZoneSize-($headerSizing-$headerSpacing))/$headerSpacing);
    my $elementZoneSize_binsDistance=($elementZoneSize_bins * $headerSpacing)+($headerSizing-$headerSpacing);
    print STDERR "\telementZoneSize\t$elementZoneSize\t$elementZoneSize_bins\t$elementZoneSize_binsDistance\n" if($verbose);

    # build headers
    my $pileUp_inc2header={};
    my $label_to_offset={};
    my $labelInc=0;
    for(my $i=-$elementZoneSize_bins;$i<=$elementZoneSize_bins;$i++) {
        my $xLabel=($i*$headerSpacing)."bp";
        $xLabel="+".$xLabel if($i > 0);
        $xLabel="x_".$xLabel;
        $pileUp_inc2header->{ x }->{$labelInc}=$xLabel;
        
        my $yLabel=($i*$headerSpacing)."bp";
        $yLabel="+".$yLabel if($i > 0);
        $yLabel="y_".$yLabel;
        $pileUp_inc2header->{ y }->{$labelInc}=$yLabel;
        
        $label_to_offset->{$labelInc}=$i;
        
        $labelInc++;
    }
    
    my $numYLabels=keys %{$pileUp_inc2header->{ y }};
    my $numXLabels=keys %{$pileUp_inc2header->{ x }};
    
    my $pileUpMatrix={};
    my $highlightMatrix={};
    
    my $numAnchors=@{$elements};
    my $numInteractions=(($numAnchors*$numAnchors)-$numAnchors);
    $numInteractions=keys %{$interactions} if( keys %{$interactions} != 0);
    my $interactionCounter=0;
    my $initFlag=1;
    my $pcComplete=0;
    
    for(my $a1=0;$a1<$numAnchors;$a1++) {
        my $element_1=$elements->[$a1]->{ name };
        my $element_1_name=$elements->[$a1]->{ name2 };
        
        my $elementObject_1=getHeaderObject($element_1,1);
        my $element1_chromosome=$elementObject_1->{ chromosome };
        
        my $elementIndex_1=-1;
        $elementIndex_1=$header2inc->{ xy }->{$element_1} if(exists($header2inc->{ xy }->{$element_1}));
        croak "non-existant header ($element_1)" if($elementIndex_1 == -1);
    
        my $elementIndex_y=$elementIndex_1;
        
        for(my $a2=0;$a2<$numAnchors;$a2++) {
        
            # exclude self:self interactions
            next if($a1 == $a2);
            
            
            my $element_2=$elements->[$a2]->{ name };
            my $element_2_name=$elements->[$a2]->{ name2 };
            my $elementObject_2=getHeaderObject($element_2,1);
            my $element2_chromosome=$elementObject_2->{ chromosome };
        
            my $key=$element_1_name."___".$element_2_name;
            next if(!exists($interactions->{$key}));
            
            my $elementIndex_2=-1;
            $elementIndex_2=$header2inc->{ xy }->{$element_2} if(exists($header2inc->{ xy }->{$element_2}));
            croak "non-existant header [$element_2]" if($elementIndex_2 == -1);
            
            my $elementIndex_x=$elementIndex_2;
            
            $interactionCounter++;
            
            my $elementDistance=abs($a1-$a2);
            my $interactionDistance=(abs($elementIndex_y-$elementIndex_x)*$headerSpacing)+($headerSizing-$headerSpacing);
            
            # subset by interaction (bp) distance
            next if((defined($minDistance)) and ($interactionDistance < $minDistance));
            next if((defined($maxDistance)) and ($interactionDistance > $maxDistance));
            
            # subset by element (#element) distance
            next if((defined($minElementDistance)) and ($elementDistance < $minElementDistance));
            next if((defined($maxElementDistance)) and ($elementDistance > $maxElementDistance));
        
            my $tmpMatrix={} if($debugMode);
            my $tmp_inc2header={} if($debugMode);
            
            # exclude diagonal 
            next if(($includeDiagonal == 0) and (abs($elementIndex_y-$elementIndex_x) <= (($elementZoneSize_bins*2)+1)));
                    
            for(my $y=0;$y<$numYLabels;$y++) {
                my $yOffset=$label_to_offset->{$y};
                my $tmpY=($elementIndex_y+$yOffset);
                
                # catch out of bounds
                next if($tmpY >= $numTotalHeaders);
                next if($tmpY < 0);
                    
                my $tmp_yHeader=$inc2header->{ xy }->{$tmpY};
                my $tmp_yHeader_object=getHeaderObject($tmp_yHeader,1);
                my $tmp_yHeader_chromosome=$tmp_yHeader_object->{ chromosome };
                
                next if($tmp_yHeader_chromosome ne $element1_chromosome);
                    
                $tmp_inc2header->{ y }->{ $y }=$tmp_yHeader if($debugMode);
                
                for(my $x=0;$x<$numXLabels;$x++) {
                    my $xOffset=$label_to_offset->{$x};
                    my $tmpX=($elementIndex_x+$xOffset);
                    
                    # catch out of bounds
                    next if($tmpX >= $numTotalHeaders);
                    next if($tmpX < 0);
                    
                    my $tmp_xHeader=$inc2header->{ xy }->{$tmpX};
                    my $tmp_xHeader_object=getHeaderObject($tmp_xHeader,1);
                    my $tmp_xHeader_chromosome=$tmp_xHeader_object->{ chromosome };
                    
                    next if($tmp_xHeader_chromosome ne $element2_chromosome);
                    
                    $tmp_inc2header->{ x }->{ $x }=$tmp_xHeader if($debugMode);
                    
                    next if($tmpX < 0);
                    
                    my $tmp_xHeaderObject=getHeaderObject($tmp_xHeader,1);
                    my $tmp_yHeaderObject=getHeaderObject($tmp_yHeader,1);
                    my $interactionDistance=getInteractionDistance($matrixObject,$tmp_yHeaderObject,$tmp_xHeaderObject,1);
                    my $interactionType=classifyInteractionDistance($interactionDistance);
                    
                    # subset by distance
                    next if((defined($minDistance)) and ($interactionDistance < $minDistance));
                    next if((defined($maxDistance)) and ($interactionDistance > $maxDistance));
                    
                    my $inten=$matrixObject->{ missingValue };
                    $inten=$matrix->{$tmpY}->{$tmpX} if(exists($matrix->{$tmpY}->{$tmpX}));
                    $inten = "NA" if($inten eq "NA");
                    
                    $highlightMatrix->{$tmpY}->{$tmpX}+=1 if($inten ne "NA");
                    $highlightMatrix->{$tmpX}->{$tmpY}+=1 if($inten ne "NA");
                    $highlightMatrix->{$tmpY}->{$tmpX}="NA" if($inten eq "NA");
                    $highlightMatrix->{$tmpX}->{$tmpY}="NA" if($inten eq "NA");
                    
                    next if($inten eq "NA");
                    
                    if(($aggregrateMode ne "sum") and ($aggregrateMode ne "mean")) {
                        @{$pileUpMatrix->{$interactionType}->{$y}->{$x}}=() if(!exists($pileUpMatrix->{$interactionType}->{$y}->{$x}));
                        push(@{$pileUpMatrix->{$interactionType}->{$y}->{$x}},$inten);
                        croak "too many data points, only sum/mean allowed [!".$aggregrateMode."]" if(@{$pileUpMatrix->{$interactionType}->{$y}->{$x}} > 10000);
                    } else {
                        $pileUpMatrix->{$interactionType}->{$y}->{$x}->{ sum }+=$inten;
                        $pileUpMatrix->{$interactionType}->{$y}->{$x}->{ count }++;
                    }
                    
                    $tmpMatrix->{$y}->{$x}=$inten if($debugMode);
                    
                }
            }
            
            $pcComplete = 100 if($interactionCounter == ($numInteractions));
            print STDERR "\e[A" if(($initFlag == 0) and ($verbose));
            printf STDERR "\t%.2f%% complete ($interactionCounter/".($numInteractions).")...\n", $pcComplete if($verbose);
            $pcComplete = round((($interactionCounter/($numInteractions))*100),2);
            $initFlag=0;
            
            my $element_1_shortName=(split(/__/,$element_1_name))[0];
            my $element_2_shortName=(split(/__/,$element_2_name))[0];
            
            my $tmpMatrixFile=$element_1_shortName."__".$element_2_shortName."__".$interactionCounter.".matrix.gz";            
            writeMatrix($tmpMatrix,$tmp_inc2header,$tmpMatrixFile,0) if($debugMode);
            
        }
    }
    
    print STDERR "\n" if($verbose);
    
    my $distanceDescriptor="";
    $distanceDescriptor="__".$minDistance."-".$maxDistance if((defined($minDistance)) or (defined($maxDistance)));
    $output .= "___".$elementFileName.$distanceDescriptor;
    
    my $highlightMatrixFile=$output.".highlight.matrix.gz";
    writeMatrix($highlightMatrix,$inc2header,$highlightMatrixFile,"NA",$commentLine);
    
    # write cis pileup matrix
    my $cis_pileUpMatrixFile=$output.".cis.pileUp.matrix.gz";
    writePileUpMatrixFile("cis",$cis_pileUpMatrixFile,$aggregrateMode,$pileUp_inc2header,$pileUpMatrix,$commentLine,$verbose) if(exists($pileUpMatrix->{ cis }) and ($excludeCis == 0));
    
    # write trans pileup matrix
    my $trans_pileUpMatrixFile=$output.".trans.pileUp.matrix.gz";
    writePileUpMatrixFile("trans",$trans_pileUpMatrixFile,$aggregrateMode,$pileUp_inc2header,$pileUpMatrix,$commentLine,$verbose) if(exists($pileUpMatrix->{ trans }) and ($excludeTrans == 0));
    
    return($cis_pileUpMatrixFile,$trans_pileUpMatrixFile);
    
}

sub writePileUpMatrixFile($$$$$$$) {
    my $interactionType=shift;
    my $outputFile=shift;
    my $aggregrateMode=shift;
    my $pileUp_inc2header=shift;
    my $numYLabels=keys %{$pileUp_inc2header->{ y }};
    my $numXLabels=keys %{$pileUp_inc2header->{ x }};
    my $pileUpMatrix=shift;
    my $commentLine=shift;
    my $verbose=shift;
    
    print STDERR "writing ($interactionType) pileUpMatrix ...\n" if($verbose);
    
    open(OUT,outputWrapper($outputFile,$commentLine)) or croak "Could not open file [$outputFile] - $!";
    
    for(my $x=0;$x<$numXLabels;$x++) {
        my $xLabel=$pileUp_inc2header->{ x }->{$x};
        print OUT "\t$xLabel";
    }
    
    print OUT "\n";
    
    for(my $y=0;$y<$numYLabels;$y++) {
        my $yLabel=$pileUp_inc2header->{ y }->{$y};
        print OUT "$yLabel";
        
        for(my $x=0;$x<$numXLabels;$x++) {
            
            my $score="NA";
            if(($aggregrateMode ne "sum") and ($aggregrateMode ne "mean")) {
                my @tmpArr=();
                @tmpArr=@{$pileUpMatrix->{ $interactionType }->{$y}->{$x}} if(exists($pileUpMatrix->{ $interactionType }->{$y}->{$x}));
                my $tmpArrSize=scalar @tmpArr;
                my $tmpArrStats=listStats(\@tmpArr) if(@tmpArr > 0);
                croak "invalid aggregrateMode [$aggregrateMode]." if(!exists($tmpArrStats->{ $aggregrateMode }));
                $score=$tmpArrStats->{ $aggregrateMode } if(exists($tmpArrStats->{ $aggregrateMode }));
            } else {
                my $sum="NA";
                $sum=$pileUpMatrix->{ $interactionType }->{$y}->{$x}->{ sum } if(exists($pileUpMatrix->{ $interactionType }->{$y}->{$x}->{ sum }));
                my $count="NA";
                $count=$pileUpMatrix->{ $interactionType }->{$y}->{$x}->{ count } if(exists($pileUpMatrix->{ $interactionType }->{$y}->{$x}->{ count }));
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
    
    print STDERR "\n" if($verbose);
}

sub buildElements($;$) {
    my $interactionBedFile=shift;
    # optional
    my $bedName="";
    $bedName=shift if @_;
    
    $bedName = "__".$bedName if($bedName ne "");
    my $combinedBedFile=$bedName."__".getSmallUniqueString().".bed";
    
    my %interactions=();
    
    open(OUT,outputWrapper($combinedBedFile)) or confess "Could not open file [$combinedBedFile] - $!";
    
    open(IN,inputWrapper($interactionBedFile)) or confess "Could not open file [$interactionBedFile] - $!";
    while(my $line = <IN>) {
        chomp($line);
        next if(($line eq "") or ($line =~ m/^#/));
        
        # skip possible BED headers
        next if(($line eq "") or ($line =~ m/^#/));
        next if($line =~ /^#/);
        next if($line =~ /^track/);
        next if($line =~ /^chrom/);
        
        my @tmp=split(/\t/,$line);
        
        confess "bad itx format - expected 7 columns, bed3 + bed3 + name\n\t@tmp\n\t$line\n" if(@tmp != 7);
        
        $tmp[1] =~ s/,//g;
        $tmp[2] =~ s/,//g;
        $tmp[4] =~ s/,//g;
        $tmp[5] =~ s/,//g;
        
        confess "ERROR-1: invalid BED format [$line]" if($line eq "");
        confess "ERROR-2: invalid BED format [$line]" if($tmp[0] !~ /^chr/);
        confess "ERROR-2: invalid BED format [$line]" if($tmp[3] !~ /^chr/);
        
        confess "ERROR-3: invalid BED format [$line]" if($tmp[1] !~ (/^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/));
        confess "ERROR-3: invalid BED format [$line]" if($tmp[4] !~ (/^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/));
        
        confess "ERROR-4: invalid BED format [$line]" if($tmp[2] !~ (/^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/));
        confess "ERROR-4: invalid BED format [$line]" if($tmp[5] !~ (/^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/));
        
        my $midpoint_1=(($tmp[1]+$tmp[2])/2);
        my $midpoint_2=(($tmp[4]+$tmp[5])/2);
        
        $tmp[1]=floor($midpoint_1);
        $tmp[2]=ceil($midpoint_1);
        
        $tmp[4]=floor($midpoint_2);
        $tmp[5]=ceil($midpoint_2);
        
        my $name_1 = $tmp[6]."_1";
        my $name_2 = $tmp[6]."_2";
        
        my $key=$name_1."___".$name_2;
        $interactions{$key}=1;
        
        print OUT "$tmp[0]\t$tmp[1]\t$tmp[2]\t$name_1\n";
        print OUT "$tmp[3]\t$tmp[4]\t$tmp[5]\t$name_2\n";
    }
    
    close(IN);
    close(OUT);
    
    return($combinedBedFile,\%interactions);
}
        
my %options;
my $results = GetOptions( \%options,'inputMatrix|i=s','verbose|v','output|o=s','elementBedFiles|ebf=s@','interactionBedFile|ibf=s','elementName|en=s','elementZoneSize|ezs=i','minDistance|minDist=i','maxDistance|maxDist=i','minElementDistance|minED=i','maxElementDistance|maxED=i','aggregrateMode|am=s','excludeZero|ez','includeDiagonal|id','excludeCis|ec','excludeTrans|et','debugMode|d') or croak help();
my ($ret,$inputMatrix,$verbose,$output,$elementBedFiles,$interactionBedFile,$elementName,$elementZoneSize,$minDistance,$maxDistance,$minElementDistance,$maxElementDistance,$aggregrateMode,$excludeZero,$includeDiagonal,$excludeCis,$excludeTrans,$debugMode)=check_options( \%options );

intro() if($verbose);

#get the absolute path of the script being used
my $cwd = getcwd();
my $fullScriptPath=abs_path($0);
my @fullScriptPathArr=split(/\//,$fullScriptPath);
@fullScriptPathArr=@fullScriptPathArr[0..@fullScriptPathArr-3];
my $scriptPath=join("/",@fullScriptPathArr);
my $commentLine=getScriptOpts($ret,$tool);

croak "inputMatrix [$inputMatrix] does not exist" if(!(-e $inputMatrix));

croak "must supply valid elementBedFile or interactionBedFile!" if((@{$elementBedFiles} == 0) and !(-e $interactionBedFile));
croak "must supply valid elementBedFile or interactionBedFile!" if((@{$elementBedFiles} != 0) and (-e $interactionBedFile));

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

my $equalHeaderFlag=$matrixObject->{ equalHeaderFlag };
croak "matrix must be symmetrical / equal sized intervals" if(($symmetrical == 0) or ($equalHeaderFlag == 0));

for(my $i=0;$i<@{$elementBedFiles};$i++) {
    my $elementBedFile = $elementBedFiles->[$i];
    print STDERR "validating $elementBedFile ...\n" if($verbose);
    validateBed($elementBedFile);
}

print STDERR "\n" if($verbose);

my $elementBedFile="";
my $elementFileName="";
my $interactions={};

if(@{$elementBedFiles} > 1) {
 
    print STDERR "using element file\n" if($verbose);
    print STDERR "\n" if($verbose);
    
    $elementBedFile=$elementBedFiles->[0] if(@{$elementBedFiles});
    $elementFileName=getFileName($elementBedFile);
    $elementFileName=$elementName if($elementName ne "");
    print STDERR "combining bed files ...\n" if($verbose);
    $elementBedFile=combineBedFiles($elementBedFiles,$elementName);
    print STDERR "\tdone\n" if($verbose);
    print STDERR "\n" if($verbose);

    # set all elements to midpoint of element
    print STDERR "translating interval bed to midpoint bed ...\n" if($verbose);
    $elementBedFile=midpointBedFile($elementBedFile,$elementName);
    print STDERR "\tdone\n" if($verbose);

    print STDERR "\n" if($verbose);

} else {

    print STDERR "building elements from itx file\n" if($verbose);
    print STDERR "\n" if($verbose);
    
    $elementFileName=getFileName($interactionBedFile);
    $elementFileName=$elementName if($elementName ne "");
    
    ($elementBedFile,$interactions)=buildElements($interactionBedFile,$elementName);
    
}

print STDERR "running headers2bed ...\n" if($verbose);
my $headerBedFile=headers2bed($matrixObject);
print STDERR "\t$headerBedFile\n" if($verbose);

print STDERR "\n" if($verbose);

print STDERR "intersecting Bed files ...\n" if($verbose);
my $bedOverlapFile=intersectBED($headerBedFile,$elementBedFile);
print STDERR "\t$bedOverlapFile\n" if($verbose);
system("rm '".$headerBedFile."'") if(-e $headerBedFile);
system("rm '".$elementBedFile."'") if(-e $elementBedFile);

print STDERR "\n" if($verbose);

print STDERR "loading Bed file ...\n" if($verbose);
my ($elements)=loadBED($bedOverlapFile);
print STDERR "\tfound ".@{$elements}." elements\n" if($verbose);
system("rm '".$bedOverlapFile."'");

croak "found no overlapping headers!" if(@{$elements} == 0);

print STDERR "\n" if($verbose);

(my $matrix,$matrixObject)=getData($inputMatrix,$matrixObject,$verbose,$minDistance,$maxDistance,$excludeCis,$excludeTrans);

print STDERR "\n" if($verbose);

# calculate the boundaryReach index for each bin and store in a new data struct.
print STDERR "piling up data ...\n" if($verbose);
interactionPileUp($matrixObject,$matrix,$elementFileName,$elements,$interactions,$elementZoneSize,$minDistance,$maxDistance,$minElementDistance,$maxElementDistance,$aggregrateMode,$excludeZero,$includeDiagonal,$excludeCis,$excludeTrans,$commentLine,$debugMode,$verbose);