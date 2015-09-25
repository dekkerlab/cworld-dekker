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

    my ($inputMatrix,$verbose,$output,$excludeCis,$excludeTrans,$skipNA,$excludeZero);
    
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
    
    croak "cannot exclude both cis and trans - no data remaining!" if(($excludeCis) and ($excludeTrans));
   
   $ret->{ inputMatrix }=$inputMatrix;
    $ret->{ verbose }=$verbose;
    $ret->{ output }=$output;
    $ret->{ excludeCis }=$excludeCis;
    $ret->{ excludeTrans }=$excludeTrans;
    $ret->{ skipNA }=$skipNA;
    $ret->{ excludeZero }=$excludeZero;
    
    return($ret,$inputMatrix,$verbose,$output,$excludeCis,$excludeTrans,$skipNA,$excludeZero);
}
 
sub intro() {
    print STDERR "\n";
    
    print STDERR "Tool:\t\t".$tool."\n";
    print STDERR "Version:\t".$cworld::dekker::VERSION."\n";
    print STDERR "Summary:\ttransform tsv matrix into 3 column tsv file\n";
    
    print STDERR "\n";
}

sub help() {
    intro();
    
    print STDERR "Usage: perl matrix2pairwise.pl [OPTIONS] -i <inputMatrix>\n";
    
    print STDERR "\n";
    
    print STDERR "Required:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-i", "[]", "input matrix file");
    
    print STDERR "\n";
    
    
    print STDERR "Options:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-v", "[]", "FLAG, verbose mode");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-o", "[]", "prefix for output file(s)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ec", "[]", "FLAG, exclude CIS data");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--et", "[]", "FLAG, exclude TRANS data");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--sna", "[]", "FLAG, skip NAs in output");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ez", "[]", "FLAG, ignore zeros from all calculations");
    
    print STDERR "\n";
    
    print STDERR "Notes:";
    print STDERR "
    This script calculates a lowess smoothing of the data per distance.
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
my $results = GetOptions( \%options,'inputMatrix|i=s','verbose|v','output|o=s','excludeCis|ec','excludeTrans|et','skipNA|sna','excludeZero|ez') or croak help();
my ($ret,$inputMatrix,$verbose,$output,$excludeCis,$excludeTrans,$skipNA,$excludeZero)=check_options( \%options );

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
my $equalHeaderFlag=$matrixObject->{ equalHeaderFlag };
my $inputMatrixName=$matrixObject->{ inputMatrixName };
$output=$matrixObject->{ output };

print STDERR "printing pairwise file...\n" if($verbose);
my $pairwiseFile=matrix2pairwise($matrixObject,$inputMatrix,$excludeCis,$excludeTrans,$skipNA,$excludeZero);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);



