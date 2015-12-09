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

    my ($inputMatrix,$verbose,$output,$minDistance,$maxDistance,$logTransform,$skipNA,$excludeZero);
    
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
    
    if( exists($opts->{ logTransform }) ) {
        $logTransform = 10;
    } else {
        $logTransform = 0;
    }
    
    if( exists($opts->{ skipNA }) ) {
        $skipNA = 1;
    } else {
        $skipNA = 0;
    }
    
    if( exists($opts->{ excludeZero }) ) {
        $excludeZero = 1;
    } else {
        $excludeZero = 0;
    }
  
    $ret->{ inputMatrix }=$inputMatrix;
    $ret->{ verbose }=$verbose;
    $ret->{ output }=$output;
    $ret->{ minDistance }=$minDistance;
    $ret->{ maxDistance }=$maxDistance;
    $ret->{ logTransform }=$logTransform;
    $ret->{ skipNA }=$skipNA;
    $ret->{ excludeZero }=$excludeZero;
    
    return($ret,$inputMatrix,$verbose,$output,$minDistance,$maxDistance,$logTransform,$skipNA,$excludeZero);
    
}

sub calculateCumulativeReadsPerDistance($$$) {
    my $cisFile=shift;
    my $matrixSum=shift;
    my $distanceFile=shift;
    
    open(OUT,outputWrapper($distanceFile)) or croak "Could not open file [$distanceFile] - $!";
    print OUT "interactionDistance\tcumulativeReads\tcumulativePercent\n";
    
    my $lastInteractionDistance=-1;
    my $cumulativeReads=0;
    
    open(IN,inputWrapper($cisFile)) or croak "Could not open file [$cisFile] - $!";
    while(my $line = <IN>) {
        chomp($line);
        next if(($line eq "") or ($line =~ m/^#/));
        
        my @tmp=split(/\t/,$line);
        
        my $interactionDistance=$tmp[0];
        my $cScore=$tmp[1];
        my $interactionKey=$tmp[2];
        
        if(($lastInteractionDistance != -1) and ($interactionDistance != $lastInteractionDistance)) {
            my $cumulativePercent=round(($cumulativeReads/$matrixSum)*100,3);
            print OUT "$lastInteractionDistance\t$cumulativeReads\t$cumulativePercent\n";
        }
        
        $cumulativeReads += $cScore;
        $lastInteractionDistance=$interactionDistance;
        
    }
    close(IN);
    
    close(OUT);
    
}

sub intro() {
    print STDERR "\n";
    
    print STDERR "Tool:\t\t".$tool."\n";
    print STDERR "Version:\t".$cworld::dekker::VERSION."\n";
    print STDERR "Summary:\tcumlative reads versus distance\n";
    
    print STDERR "\n";
}

sub help() {
    intro();
    
    print STDERR "Usage: perl matrix2distance.pl [OPTIONS] -i <inputMatrix>\n";
    
    print STDERR "\n";
    
    print STDERR "Required:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-i", "[]", "input matrix file");
    
    print STDERR "\n";
    
    print STDERR "Options:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-v", "[]", "FLAG, verbose mode");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-o", "[]", "prefix for output file(s)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--minDist", "[]", "minimum allowed interaction distance, exclude < N distance (in BP)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--maxDist", "[]", "maximum allowed interaction distance, exclude > N distance (in BP)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--lt", "[]", "optional, log transform data (2,10)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ez", "[]", "FLAG, exclude zeros from all calculations");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--sn", "[]", "FLAG, skip NAs in calculation");
    print STDERR "\n";
    
    print STDERR "Notes:";
    print STDERR "
    This script can compute the cumlative reads across distances for a given matrix.
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
my $results = GetOptions( \%options,'inputMatrix|i=s','verbose|v','output|o=s','minDistance|minDist=i','maxDistance|maxDist=i','logTransform|lt','skipNA|sn','excludeZero|ez') or croak help();
my ($ret,$inputMatrix,$verbose,$output,$minDistance,$maxDistance,$logTransform,$skipNA,$excludeZero)=check_options( \%options );

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

my $excludeCis=0;
my $excludeTrans=1;

print STDERR "printing pairwise-distance file...\n" if($verbose);
my $pairwiseFile=matrix2distance($matrixObject,$inputMatrix,$excludeCis,$excludeTrans,$skipNA,$excludeZero);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

print STDERR "seperating cis/trans data...\n" if($verbose);
my ($cisFile,$transFile,$matrixSum)=matrix2listfile($matrixObject,$inputMatrix,$verbose,$excludeZero,$minDistance,$maxDistance,$excludeCis,$excludeTrans);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

print STDERR "calculating cumlative read distribution...\n" if($verbose);
my $distanceFile=$output.".cum.distance.txt.gz";
calculateCumulativeReadsPerDistance($cisFile,$matrixSum,$distanceFile);
print STDERR "\tdone\n" if($verbose);

system("rm '".$cisFile."'") if(-e($cisFile));
system("rm '".$transFile."'") if(-e($transFile));

print STDERR "\n" if($verbose);

# print loess
print STDERR "plotting distance distribution...\n" if($verbose);
system("Rscript '".$scriptPath."/R/plotDistance.R' '".$cwd."' '".$distanceFile."' ".$output." ".$logTransform." > /dev/null");
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);
