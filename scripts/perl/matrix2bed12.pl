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

    my ($inputMatrix,$verbose,$output,$scaleStart,$scaleEnd,$scaleStartTile,$scaleEndTile);
    
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
    
     if( exists($opts->{ scaleStart }) ) {
        $scaleStart = $opts->{ scaleStart };
    } else {
        $scaleStart = "NA";
    }
    
    if( exists($opts->{ scaleEnd }) ) {
        $scaleEnd = $opts->{ scaleEnd };
    } else {
        $scaleEnd = "NA";
    }
    
    if( exists($opts->{ scaleStartTile }) ) {
        $scaleStartTile = $opts->{ scaleStartTile };
    } else {
        $scaleStartTile = 0.025;
    }
    
    if( exists($opts->{ scaleEndTile }) ) {
        $scaleEndTile = $opts->{ scaleEndTile };
    } else {
        $scaleEndTile = 0.975;
    }
    
    $ret->{ inputMatrix }=$inputMatrix;
    $ret->{ verbose }=$verbose;
    $ret->{ output }=$output;
    $ret->{ scaleStart }=$scaleStart;
    $ret->{ scaleEnd }=$scaleEnd;
    $ret->{ scaleStartTile }=$scaleStartTile;
    $ret->{ scaleEndTile }=$scaleEndTile;
    
    return($ret,$inputMatrix,$verbose,$output,$scaleStart,$scaleEnd,$scaleStartTile,$scaleEndTile);
}

sub intro() {
    print STDERR "\n";
    
    print STDERR "Tool:\t\t".$tool."\n";
    print STDERR "Version:\t".$cworld::dekker::VERSION."\n";
    print STDERR "Summary:\ttransform matrix into bed12 format (track per row)\n";
    
    print STDERR "\n";
}

sub help() {
    intro();
    
    print STDERR "Usage: perl matrix2bed12.pl [OPTIONS] -i <inputMatrix>\n";
    
    print STDERR "\n";
    
    print STDERR "Required:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-i", "[]", "input matrix file");
    
    print STDERR "\n";
    print STDERR "Options:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-v", "[]", "FLAG, verbose mode");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-o", "[]", "prefix for output file(s)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--start", "[]", "absolute value for display start");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--end", "[]", "absolute vlaue for display end");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--startTile", "[0.025]", "fraction value for display start");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--endTile", "[0.975]", "fraction value for display end");
    
    print STDERR "\n";
    
    print STDERR "Notes:";
    print STDERR "
    This script transform a matrix into bed12 format (1 track per row/col).
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

sub matrix2bed12($$$$$$) {
    my $matrixObject=shift;
    my $matrix=shift;
    my $scaleStart=shift;
    my $scaleEnd=shift;
    my $scaleBucketSize=shift;
    my $nColorShades=shift;
    
    my $inc2header=$matrixObject->{ inc2header };
    my $numYHeaders=$matrixObject->{ numYHeaders };
    my $numXHeaders=$matrixObject->{ numXHeaders };
    my $missingValue=$matrixObject->{ missingValue };
    my $symmetrical=$matrixObject->{ symmetrical };
    my $output=$matrixObject->{ output };
    
    my %usedName=();
    
    my $bedFile=$output.".bed.gz";
    open(OUT,outputWrapper($bedFile)) or croak "Could not open file [$bedFile] - $!";
    print OUT "track name=$output description='$output ($scaleStart - $scaleEnd)' useScore=1 visibility=4 itemRgb='On'\n";

    for(my $y=0;$y<$numYHeaders;$y++) {
        my $yHeader=$inc2header->{ y }->{$y};
        my $yHeaderObject=getHeaderObject($yHeader);
        my $yChromosome=$yHeaderObject->{ chromosome };
        for(my $x=0;$x<$numXHeaders;$x++) {
            my $xHeader=$inc2header->{ x }->{$x};    
            my $xHeaderObject=getHeaderObject($xHeader);
            my $xChromosome=$yHeaderObject->{ chromosome };
            
            my $interactionDistance=getInteractionDistance($matrixObject,$yHeaderObject,$xHeaderObject,1);
            
            #cis only
            next if($interactionDistance == -1);
            next if($yChromosome ne $xChromosome);
            
            my $score=$matrixObject->{ missingValue };
            $score=$matrix->{$y}->{$x} if(exists($matrix->{$y}->{$x}));
            
            next if($score == 0);
            next if($score eq "NA");
            
            my $scaleIndex = -1;
            $scaleIndex=getColorIndex($score,$scaleStart,$scaleEnd,$nColorShades,$scaleBucketSize) if($score ne "NA");
            next if($scaleIndex == 0);
            
            my $chromosome=$yChromosome=$xChromosome;
            
            my $yStart=$yHeaderObject->{ start };
            my $yEnd=$yHeaderObject->{ end };
            my $xStart=$xHeaderObject->{ start };
            my $xEnd=$xHeaderObject->{ end };
            
            my $start=min($yStart,$xStart);
            my $end=max($yEnd,$xEnd);
            
            my $name=$yHeader.".".$xHeader;
            my $name2=$xHeader.".".$yHeader;
            
            next if((exists($usedName{$name})) or (exists($usedName{$name2})));
            $usedName{$name}=1;
            $usedName{$name2}=1;
            
            my $strand=".";
            my $rgb="0,0,0";
            $rgb="255,0,0" if($score > 0);
            $rgb="0,0,255" if($score < 0);
            
            my $blockCount=2;
            
            my ($blockSizes,$blockStarts);
            if($yStart < $xStart) {
                $blockSizes=($yEnd-$yStart).",".($xEnd-$xStart);
                $blockStarts="0,".($xStart-$start);
            } else {
                $blockSizes=($xEnd-$xStart).",".($yEnd-$yStart);
                $blockStarts="0,".($yStart-$start);
            }
        
            print OUT "$chromosome\t$start\t$end\t$name\t$scaleIndex\t$strand\t$start\t$end\t$rgb\t$blockCount\t$blockSizes\t$blockStarts\n";
            
        }
    }

    close(OUT);
    
    return($bedFile);
}

my %options;
my $results = GetOptions( \%options,'inputMatrix|i=s','verbose|v','output|o=s','scaleStart|start=s','scaleEnd|end=s','scaleStartTile|startTile=s','scaleEndTile|endTile=s') or croak help();
my ($ret,$inputMatrix,$verbose,$output,$scaleStart,$scaleEnd,$scaleStartTile,$scaleEndTile)=check_options( \%options );

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

#read Matrix
my ($matrix)=getData($inputMatrix,$matrixObject,$verbose);

print STDERR "\n" if($verbose);

my $nColorShades=1000;
print STDERR "Calculating Auto Color scale [$scaleStartTile:$scaleEndTile]...\n" if($verbose);
my ($tmpScaleStart,$tmpScaleEnd)=autoScale($matrixObject,$matrix,$scaleStartTile,$scaleEndTile);
$scaleStart=$tmpScaleStart if($scaleStart eq "NA");
$scaleEnd=$tmpScaleEnd if($scaleEnd eq "NA");
my $scaleDistance = ($scaleEnd-$scaleStart);
my $scaleBucketSize = ($scaleDistance/($nColorShades-1));
print STDERR "\tALL\t".$scaleStart."-".$scaleEnd."\n" if($verbose);

print STDERR "\n" if($verbose);

print STDERR "Transforming into bed12 ...\n" if($verbose);
my $bedFile=matrix2bed12($matrixObject,$matrix,$scaleStart,$scaleEnd,$scaleBucketSize,$nColorShades);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);