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
    
    my ($inputMatrix,$verbose,$output,$xHeaderFile,$yHeaderFile);
    
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
    
    if( exists($opts->{ xHeaderFile }) ) {
        $xHeaderFile = $opts->{ xHeaderFile };
    } else {
        print STDERR "\nERROR: Option xHeaderFile|xhf is required.\n";
        help();
    }
    
    if( exists($opts->{ yHeaderFile }) ) {
        $yHeaderFile = $opts->{ yHeaderFile };
    } else {
        print STDERR "\nERROR: Option yHeaderFile|yhf is required.\n";
        help();
    }
    
    $ret->{ inputMatrix }=$inputMatrix;
    $ret->{ verbose }=$verbose;
    $ret->{ output }=$output;
    $ret->{ xHeaderFile }=$xHeaderFile;
    $ret->{ yHeaderFile }=$yHeaderFile;
    
    return($ret,$inputMatrix,$verbose,$output,$xHeaderFile,$yHeaderFile);
}

sub intro() {
    print STDERR "\n";
    
    print STDERR "Tool:\t\t".$tool."\n";
    print STDERR "Version:\t".$cworld::dekker::VERSION."\n";
    print STDERR "Summary:\tadd headers to a matrix txt file\n";
    
    print STDERR "\n";
}

sub help() {
    intro();
    
    print STDERR "Usage: perl addMatrixHeaders.pl [OPTIONS] -i <inputMatrix> -xhf <xHeaderFile> -yhf <yHeaderFile>\n";
    
    print STDERR "\n";
    
    print STDERR "Required:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-i", "[]", "input matrix file");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--xhf", "[]", "x-axis header file");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--yhf", "[]", "y-axis header file");
    
    
    print STDERR "\n";
    
    print STDERR "Options:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-v", "[]", "FLAG, verbose mode");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-o", "[]", "prefix for output file(s)");
    
    print STDERR "\n";
    
    print STDERR "Notes:";
    print STDERR "
    This script can add headers to a _naked_ matrix file.
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

sub loadHeaders($) {
    my $headerFile=shift;
    
    croak "headerFile [$headerFile] does not exist" if(!(-e $headerFile));
    
    my @headerArr=();
    open(IN,inputWrapper($headerFile)) or croak "Could not open file [$headerFile] - $!";
    while(my $line = <IN>) {
        chomp($line);
        next if(($line eq "") or ($line =~ m/^#/));
        
        push(@headerArr,$line);
        
    }
    close(IN);
    
    return(@headerArr);
}

sub addMatrixHeaders($$$$$$) {
    my $matrixObject=shift;
    my $inputMatrix=shift;
    my $outputPrefix=shift;
    my $xHeaderFile=shift;
    my $yHeaderFile=shift;
    my $commentLine=shift;
    
    my $headeredMatrix=$outputPrefix.".addedHeaders.matrix.gz";
    open(OUT,outputWrapper($headeredMatrix,$commentLine)) or croak "Could not open file [$headeredMatrix] - $!";
    
    my @yHeaderArr=loadHeaders($yHeaderFile);
    my @xHeaderArr=loadHeaders($xHeaderFile);
    
    my $lineNum=0;
    open(IN,inputWrapper($inputMatrix)) or croak "Could not open file [$inputMatrix] - $!";
    while(my $line = <IN>) {
        chomp($line);
        next if(($line eq "") or ($line =~ m/^#/));
        
        if($lineNum == 0) {
            print OUT "\t".join("\t",@xHeaderArr)."\n";
        } 
        print OUT $yHeaderArr[$lineNum]."\t".$line."\n";
        $lineNum++;
    }
    close(IN);
    
    close(OUT);
    
}

my %options;
my $results = GetOptions( \%options,'inputMatrix|i=s','verbose|v','output|o=s','xHeaderFile|xhf=s','yHeaderFile|yhf=s') or croak help();
my ($ret,$inputMatrix,$verbose,$output,$xHeaderFile,$yHeaderFile)=check_options( \%options );

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

print STDERR "adding matrix headers ...\n" if($verbose);
addMatrixHeaders($matrixObject,$inputMatrix,$output,$xHeaderFile,$yHeaderFile,$commentLine);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);
