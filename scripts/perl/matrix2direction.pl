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

    my ($inputMatrix,$verbose,$output,$directionMode,$directionSize,$excludeZero,$yBound,$transparentBGFlag);
    
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
    
    if( exists($opts->{ directionMode }) ) {
        $directionMode = $opts->{ directionMode };
        croak "invalid direction mode! ($directionMode) [mean,median]" if(($directionMode ne "mean") and ($directionMode ne "median"));
    } else {
        $directionMode = "mean"
    }
    
    if( exists($opts->{ directionSize }) ) {
        $directionSize = $opts->{ directionSize };
    } else {
        $directionSize = 500000;
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
    
    $ret->{ inputMatrix }=$inputMatrix;
    $ret->{ verbose }=$verbose;
    $ret->{ output }=$output;
    $ret->{ directionMode }=$directionMode;
    $ret->{ directionSize }=$directionSize;
    $ret->{ excludeZero }=$excludeZero;
    $ret->{ yBound }=$yBound;
    $ret->{ transparentBGFlag }=$transparentBGFlag;
    
    return($ret,$inputMatrix,$verbose,$output,$directionMode,$directionSize,$excludeZero,$yBound,$transparentBGFlag);
}


sub intro() {
    print STDERR "\n";
    
    print STDERR "Tool:\t\t".$tool."\n";
    print STDERR "Version:\t".$cworld::dekker::VERSION."\n";
    print STDERR "Summary:\tcalculate directionality [tads] on matrix\n";
    
    print STDERR "\n";
}

sub help() {
    intro();
    
    print STDERR "Usage: perl matrix2direction.pl [OPTIONS] -i <inputMatrix>\n";
    
    print STDERR "\n";
    
    print STDERR "Required:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-i", "[]", "input matrix file");
    
    print STDERR "\n";
    print STDERR "Options:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-v", "[]", "FLAG, verbose mode");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-o", "[]", "prefix for output file(s)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--dm", "[mean]", "optional, directionality aggregrate mode [mean,median]");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ds", "[500000]", "optional, directionality window size in bp");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ez", "[]", "FLAG, ignore 0s in all calculations");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--yb", "[auto]", "optional, -yBound - +yBound for direction plot y-axis range");    
    printf STDERR ("\t%-10s %-10s %-10s\n", "--bg", "[]", "FLAG, use transparent background of direction plot");
    print STDERR "\n";
    
    print STDERR "Notes:";
    print STDERR "
    This script can calculate the directionality log2[(eft/right) for each locus.
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

sub calculateDirection($$$$$) {
    my $matrixObject=shift;
    my $matrix=shift;
    my $directionMode=shift;
    my $directionSize=shift;
    my $directionSizeDistance=shift;
    
    my $inc2header=$matrixObject->{ inc2header };
    my $header2inc=$matrixObject->{ header2inc };
    my $numYHeaders=$matrixObject->{ numYHeaders };
    my $numXHeaders=$matrixObject->{ numXHeaders };    
    my $numHeaders=$matrixObject->{ numTotalHeaders };    
    my $symmetrical=$matrixObject->{ symmetrical };

    my %matrixDirection=();
    
    my $directionStartIndex=$directionSize;
    my $directionEndIndex=($numHeaders-$directionSize);
    
    for(my $y=$directionStartIndex;$y<$directionEndIndex;$y++) {
        my $yHead=$inc2header->{ y }->{$y};
        
        my $yObject=getHeaderObject($yHead);
        my $yChromosome=$yObject->{ chromosome };
        my $yStart=$yObject->{ start };
        my $yEnd=$yObject->{ end };
        my $yMid=$yObject->{ midpoint };
        
        my @leftData=();
        my @rightData=();
                
        # get left data
        for(my $x=0;$x<$numXHeaders;$x++) {
            my $xHead=$inc2header->{ x }->{$x};
            
            my $xObject=getHeaderObject($xHead);
            my $xChromosome=$xObject->{ chromosome };
            my $xStart=$xObject->{ start };
            my $xEnd=$xObject->{ end };
            my $xMid=$xObject->{ midpoint };
            
            # only keep those X that are 5' of Y
            next if($xMid > $yMid);
            next if(abs($xMid - $yMid) > $directionSizeDistance);
            
            #read Matrix
            my $inten=$matrixObject->{ missingValue };
            $inten=$matrix->{$y}->{$x} if(exists($matrix->{$y}->{$x}));
            next if($inten eq "NA");
            
            push(@leftData,$inten);            
        }
        
        # get right data
        for(my $x=0;$x<$numXHeaders;$x++) {
            my $xHead=$inc2header->{ x }->{$x};
            
            my $xObject=getHeaderObject($xHead);
            my $xChromosome=$xObject->{ chromosome };
            my $xStart=$xObject->{ start };
            my $xEnd=$xObject->{ end };
            my $xMid=$xObject->{ midpoint };

            # only keep those X that are 3' of Y
            next if($xMid < $yMid);
            next if(abs($xMid - $yMid) > $directionSizeDistance);
            
            my $inten=$matrixObject->{ missingValue };
            $inten=$matrix->{$y}->{$x} if(exists($matrix->{$y}->{$x}));
            next if($inten eq "NA");
            
            push(@rightData,$inten);
        }
        
        my $leftDataStats=listStats(\@leftData) if(@leftData > 0);
        my $leftDirection="NA";
        $leftDirection=$leftDataStats->{ $directionMode } if(@leftData > 0);
        
        my $rightDataStats=listStats(\@rightData) if(@rightData > 0);
        my $rightDirection="NA";
        $rightDirection=$rightDataStats->{ $directionMode } if(@rightData > 0);
        
        my $directionality="NA";
        $directionality=log($leftDirection/$rightDirection)/log(2) if(($leftDirection ne "NA") and ($rightDirection ne "NA"));        
        
        $matrixDirection{$yHead}{ direction }=$directionality;
    
    }
    
    if($symmetrical == 0) {
        
        for(my $x=0;$x<$numXHeaders;$x++) {
            my $xHead=$inc2header->{ x }->{$x};
            
            next if($x < 50);
            next if($x > ($numXHeaders-50));
            
            my $xObject=getHeaderObject($xHead);
            my $xChromosome=$xObject->{ chromosome };
            my $xStart=$xObject->{ start };
            my $xEnd=$xObject->{ end };
            my $xMid=$xObject->{ midpoint };
            
            my @leftData=();
            my @rightData=();
                        
            # get left data
            for(my $y=0;$y<$numYHeaders;$y++) {
                my $yHead=$inc2header->{ y }->{$y};
                
                my $yObject=getHeaderObject($yHead);
                my $yChromosome=$yObject->{ chromosome };
                my $yStart=$yObject->{ start };
                my $yEnd=$yObject->{ end };
                my $yMid=$yObject->{ midpoint };
                
                # only keep those Y that are 5' of X
                next if($yMid > $xMid);
                
                my $inten=$matrixObject->{ missingValue };
                $inten=$matrix->{$y}->{$x} if(exists($matrix->{$y}->{$x}));
                next if($inten eq "NA");
                
                push(@leftData,$inten);            
            }
            
            # get right data
            for(my $y=0;$y<$numYHeaders;$y++) {
                my $yHead=$inc2header->{ y }->{$y};
                
                my $yObject=getHeaderObject($yHead);
                my $yChromosome=$yObject->{ chromosome };
                my $yStart=$yObject->{ start };
                my $yEnd=$yObject->{ end };
                my $yMid=$yObject->{ midpoint };

                # only keep those Y that are 3' of X
                next if($yMid < $xMid);
                
                my $inten=$matrixObject->{ missingValue };
                $inten=$matrix->{$y}->{$x} if(exists($matrix->{$y}->{$x}));
                next if($inten eq "NA");
                
                push(@rightData,$inten);
            }
            
            my $leftDataStats=listStats(\@leftData) if(@leftData > 0);
            my $leftDirection="NA";
            $leftDirection=$leftDataStats->{ $directionMode } if(@leftData > 0);
            
            my $rightDataStats=listStats(\@rightData) if(@rightData > 0);
            my $rightDirection="NA";
            $rightDirection=$rightDataStats->{ $directionMode } if(@rightData > 0);
                    
            my $directionality="NA";
            $directionality=log($leftDirection/$rightDirection)/log(2) if(($leftDirection ne "NA") and ($rightDirection ne "NA"));
            
            $matrixDirection{$xHead}{ direction }=$directionality;
        }
    
    }
        
    return(\%matrixDirection);
    
}

sub outputDirection($$$$) {
    my $matrixObject=shift;
    my $matrixDirection=shift;
    my $directionFileName=shift;
    my $yBound=shift;
    
    my $inc2header=$matrixObject->{ inc2header };
    my $header2inc=$matrixObject->{ header2inc };
    my $numHeaders=$matrixObject->{ numYHeaders };
    my $headerSizing=$matrixObject->{ headerSizing };
    my $headerSpacing=$matrixObject->{ headerSpacing };
    
    # main direction file
    open(OUT,outputWrapper($directionFileName)) or croak "Could not open file [$directionFileName] - $!";
    print OUT "header\tstart\tend\tmidpoint\tbinStart\tbinEnd\tbinMidpoint\tdirectionScore\n";
    
    # bed graph file of direction data
    open(BEDGRAPH,outputWrapper($directionFileName.".bedGraph")) or croak "Could not open file [$directionFileName] - $!";
    print BEDGRAPH "track type=bedGraph name='".$directionFileName."' description='".$directionFileName." - insutation score' maxHeightPixels=128:64:32 visibility=full autoScale=off viewLimits=-".$yBound.":".$yBound." color=0,0,0 altColor=100,100,100\n" if($yBound != 0);
    print BEDGRAPH "track type=bedGraph name='".$directionFileName."' description='".$directionFileName." - insutation score' maxHeightPixels=128:64:32 visibility=full autoScale=on color=0,0,0 altColor=100,100,100\n" if($yBound == 0);
    
    for(my $y=0;$y<$numHeaders;$y++) {
        # dump direction data to file
        my $yHead=$inc2header->{ y }->{$y};
        
        my $direction="NA";
        $direction=$matrixDirection->{$yHead}{ direction } if(exists($matrixDirection->{$yHead}{ direction }));
        
        my $yHeadObject=getHeaderObject($yHead);
        my $yHeadChromosome=$yHeadObject->{ chromosome };
        my $yHeadStart=$yHeadObject->{ start };
        my $yHeadEnd=$yHeadObject->{ end };
        my $yHeadMidpoint=round(($yHeadStart+$yHeadEnd)/2);
        
        my $binStart = round($yHeadStart/$headerSpacing);
        my $binEnd = round($yHeadEnd/$headerSpacing);
        my $binMidpoint=(($binStart+$binEnd)/2);
        
        print OUT "$yHead\t$yHeadStart\t$yHeadEnd\t$yHeadMidpoint\t$binStart\t$binEnd\t$binMidpoint\t$direction\n";
        
        # strip off chr group if exists for proper UCSC usage
        $yHeadChromosome=stripChromosomeGroup($yHeadChromosome);
        
        $direction=0 if($direction eq "NA");
        print BEDGRAPH "$yHeadChromosome\t$yHeadStart\t$yHeadEnd\t$direction\n";
    }
    
    close(OUT);
    close(BEDGRAPH);
    
}

my %options;
my $results = GetOptions( \%options,'inputMatrix|i=s','verbose|v','output|o=s','directionMode|dm=s','directionSize|ds=i','excludeZero|ez','yBound|yb=f','transparentBGFlag|bg') or croak help();
my ($ret,$inputMatrix,$verbose,$output,$directionMode,$directionSize,$excludeZero,$yBound,$transparentBGFlag)=check_options( \%options );

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
my $contigList=$matrixObject->{ contigList };
my $numContigs=$matrixObject->{ numContigs };
my $index2contig=$matrixObject->{ index2contig };
my $symmetrical=$matrixObject->{ symmetrical };
my $inputMatrixName=$matrixObject->{ inputMatrixName };
$output=$matrixObject->{ output };

croak "matrix must be symmetrical" if($symmetrical == 0);

my $imageWidth=$numTotalHeaders*2;
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

my $numHeaders=$numYHeaders=$numXHeaders;

# get fragment spacing (i.e. bin size)
my ($equalSpacingFlag,$equalSizingFlag,$headerSpacing,$headerSizing)=getHeaderSpacing($inc2header->{ y });

my $idealBinSize=$directionSize;
$directionSize=ceil(($idealBinSize-($headerSizing-$headerSpacing))/$headerSpacing);
$directionSize = 2 if($directionSize <= 1);
my $directionSizeDistance=($directionSize * $headerSpacing)+($headerSizing-$headerSpacing);
print STDERR "directionSize\t$directionSizeDistance bp\t$directionSize\n" if($verbose);

$output .= "--dm".$directionMode;
$output .= "--ds".$directionSizeDistance;

my ($matrix)=getData($inputMatrix,$matrixObject,$verbose);

print STDERR "\n" if($verbose);

# calculate the direction index for each bin and store in a new data struct.
print STDERR "calculating direction index...\n" if($verbose);
my ($matrixDirection)=calculateDirection($matrixObject,$matrix,$directionMode,$directionSize,$directionSizeDistance);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

# write direction data to file
my $directionFileName=$output.".direction.gz";
print STDERR "outputing direction...\n" if($verbose);
outputDirection($matrixObject,$matrixDirection,$directionFileName,$yBound);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

# plot direction using R
print STDERR "plotting direction...\n" if($verbose);
system("Rscript '".$scriptPath."/R/matrix2direction.R' '".$cwd."' '".$directionFileName."' ".$headerSizing." ".$headerSpacing." ".$imageWidth." ".$yBound." ".$transparentBGFlag." > /dev/null");
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);