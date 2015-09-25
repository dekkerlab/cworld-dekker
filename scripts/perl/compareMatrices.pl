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

    my ($inputMatrix_1,$inputMatrix_2,$verbose,$output,$compareMode);
    
    my $ret={};
    
    if( exists($opts->{ inputMatrix_1 }) ) {
        $inputMatrix_1 = $opts->{ inputMatrix_1 };
    } else {
        print STDERR "\nERROR: Option inputMatrix_1|1 is required.\n";
        help();
    }
    
    if( exists($opts->{ inputMatrix_2 }) ) {
        $inputMatrix_2 = $opts->{ inputMatrix_2 };
    } else {
        print STDERR "\nERROR: Option inputMatrix_1|1 is required.\n";
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
    
    if( exists($opts->{ compareMode }) ) {
        $compareMode = $opts->{ compareMode };
    } else {
        $compareMode = "log2ratio";
    }
    
    $ret->{ inputMatrix_1 }=$inputMatrix_1;
    $ret->{ inputMatrix_2 }=$inputMatrix_2;
    $ret->{ verbose }=$verbose;
    $ret->{ output }=$output;
    $ret->{ compareMode }=$compareMode;
    
    return($ret,$inputMatrix_1,$inputMatrix_2,$verbose,$output,$compareMode);
}

sub intro() {
    print STDERR "\n";
    
    print STDERR "Tool:\t\t".$tool."\n";
    print STDERR "Version:\t".$cworld::dekker::VERSION."\n";
    print STDERR "Summary:\tperforms comparison between two matrices\n";
    
    print STDERR "\n";
}

sub help() {
    intro();
    
    print STDERR "Usage: perl compareMatrices.pl [OPTIONS] -1 <inputMatrix_1> -2 <inputMatrix_2>\n";
    
    print STDERR "\n";
    
    print STDERR "Required:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-1", "[]", "input matrix 1 file");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-2", "[]", "input matrix 2 file");
    
    print STDERR "\n";
    
    print STDERR "Options:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-v", "[]", "FLAG, verbose mode");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-o", "[]", "prefix for output file(s)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--cm", "[]", "optional, compare mode [log2ratio,add,sum,mean,subtract,divide,multiply,min,max]");

    print STDERR "\n";
    
    print STDERR "Notes:";
    print STDERR "
    This script can compare two matrices.
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
my $results = GetOptions( \%options,'inputMatrix_1|1=s','inputMatrix_2|2=s','verbose|v','output|o=s','compareMode|cm=s') or croak help();
my ($ret,$inputMatrix_1,$inputMatrix_2,$verbose,$output,$compareMode)=check_options( \%options );

intro() if($verbose);

#get the absolute path of the script being used
my $cwd = getcwd();
my $fullScriptPath=abs_path($0);
my @fullScriptPathArr=split(/\//,$fullScriptPath);
@fullScriptPathArr=@fullScriptPathArr[0..@fullScriptPathArr-3];
my $scriptPath=join("/",@fullScriptPathArr);

# ensure input matrices exist
croak "inputMatrix [$inputMatrix_1] does not exist" if(!(-e $inputMatrix_1));
croak "inputMatrix [$inputMatrix_2] does not exist" if(!(-e $inputMatrix_2));
croak "inputMatrix_1 == inputMatrix_2" if($inputMatrix_1 eq $inputMatrix_2);

my $inputMatrixName_1=getFileName($inputMatrix_1);
my $inputMatrixName_2=getFileName($inputMatrix_2);
$output=$inputMatrixName_1."___".$inputMatrixName_2 if($output eq "");

#validating headers
print STDERR "validating headers ...\n" if($verbose);
validateIdenticalMatrixStructure($inputMatrix_1,$inputMatrix_2);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

# get matrix information
my $matrixObject_1=getMatrixObject($inputMatrix_1,$output,$verbose);
my $inc2header_1=$matrixObject_1->{ inc2header };
my $header2inc_1=$matrixObject_1->{ header2inc };
my $numYHeaders_1=$matrixObject_1->{ numYHeaders };
my $numXHeaders_1=$matrixObject_1->{ numXHeaders };

my $matrix_1={};
($matrix_1)=getData($inputMatrix_1,$matrixObject_1,$verbose);

print STDERR "\n" if($verbose);

# get matrix information
my $matrixObject_2=getMatrixObject($inputMatrix_2,$output,$verbose);
my $inc2header_2=$matrixObject_2->{ inc2header };
my $header2inc_2=$matrixObject_2->{ header2inc };
my $numYHeaders_2=$matrixObject_2->{ numYHeaders };
my $numXHeaders_2=$matrixObject_2->{ numXHeaders };

my $matrix_2={};
($matrix_2)=getData($inputMatrix_2,$matrixObject_2,$verbose);

print STDERR "\n" if($verbose);

# assume both 1 and 2 are equal
my $inc2header=$inc2header_1;
my $header2inc=$header2inc_1;

# do user specified compare mode
print STDERR "calculating ($compareMode) of matrices...\n" if($verbose);
my $compareMatrix=compareMatrices($matrixObject_1,$matrixObject_2,$matrix_1,$matrix_2,$inc2header,$compareMode);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

my $compareMatrixFile=$output.".".$compareMode.".matrix.gz";
print STDERR "writing compare matrix ($compareMatrixFile)...\n" if($verbose);
writeMatrix($compareMatrix,$inc2header,$compareMatrixFile,"NA");
print STDERR "\tdone\n" if($verbose);
    
print STDERR "\n" if($verbose);
