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

    my ($inputMatrix,$verbose,$output,$headerMapFile);
    
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
    
    if( exists($opts->{ headerMapFile }) ) {
        $headerMapFile = $opts->{ headerMapFile };
    } else {
        print STDERR "\nERROR: Option headerMapFile|hmf is required.\n";
        help();
    }
    
    $ret->{ inputMatrix }=$inputMatrix;
    $ret->{ verbose }=$verbose;
    $ret->{ output }=$output;
    $ret->{ headerMapFile }=$headerMapFile;
    
    return($ret,$inputMatrix,$verbose,$output,$headerMapFile);
}

sub intro() {
    print STDERR "\n";
    
    print STDERR "Tool:\t\t".$tool."\n";
    print STDERR "Version:\t".$cworld::dekker::VERSION."\n";
    print STDERR "Summary:\treplace matrix row/col headers, or subset matrix by list of headers\n";
    
    print STDERR "\n";
}

sub help() {
    intro();
    
    print STDERR "Usage: perl changeMatrixHeaders.pl [OPTIONS] -i <inputMatrix> -hmf <headerMapFile>\n";
    
    print STDERR "\n";
    
    print STDERR "Required:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-i", "[]", "input matrix file");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--hmf", "[]", "header map file, old (tab) new");
    
    print STDERR "\n";
    
    print STDERR "Options:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-v", "[]", "FLAG, verbose mode");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-o", "[]", "prefix for output file(s)");
    
    print STDERR "\n";
    
    print STDERR "Notes:";
    print STDERR "
    This script can bin a matrix into fixed size intervals and aggregrate the data.
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

sub loadNewHeaders($$) {
    my $matrixObject=shift;
    my $headerMapFile=shift;
    
    croak "headerMapFile [$headerMapFile] does not exist!" if(!(-e $headerMapFile));
    
    my $inc2header=$matrixObject->{ inc2header };
    my $header2inc=$matrixObject->{ header2inc };
    my $numYHeaders=$matrixObject->{ numYHeaders };
    my $numXHeaders=$matrixObject->{ numXHeaders };

    open(IN,inputWrapper($headerMapFile)) or croak "Could not open file [$headerMapFile] - $!";
    while(my $line = <IN>) {
        chomp($line);
        next if(($line eq "") or ($line =~ m/^#/));
        
        my $oldHeader="NA";
        my $newHeader="NA";
        
        my @tmp=split(/\t/,$line);
        
        $oldHeader=$tmp[0] if(@tmp > 0);
        $newHeader=$oldHeader if(@tmp > 0);
        $newHeader=$tmp[1] if(@tmp > 1);
        
        if(exists($header2inc->{ x }->{$oldHeader})) {
            my $headerInc=$header2inc->{ x }->{$oldHeader};
            
            delete($header2inc->{ x }->{$oldHeader});
            delete($inc2header->{ x }->{$headerInc});
            
            $header2inc->{ x }->{$newHeader}=$headerInc;
            $inc2header->{ x }->{$headerInc}=$newHeader;
        }
        
        if(exists($header2inc->{ y }->{$oldHeader})) {
            my $headerInc=$header2inc->{ y }->{$oldHeader};
            
            delete($header2inc->{ y }->{$oldHeader});
            delete($inc2header->{ y }->{$headerInc});
            
            $header2inc->{ y }->{$newHeader}=$headerInc;
            $inc2header->{ y }->{$headerInc}=$newHeader;
        }
        
    }
    
    return($inc2header,$header2inc);
}
    
my %options;
my $results = GetOptions( \%options,'inputMatrix|i=s','verbose|v','output|o=s','headerMapFile|hmf=s') or croak help();
my ($ret,$inputMatrix,$verbose,$output,$headerMapFile)=check_options( \%options );

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

#read Matrix
my $matrix={};
($matrix)=getData($inputMatrix,$matrixObject,$verbose);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

my ($new_inc2header,$new_header2inc)=loadNewHeaders($matrixObject,$headerMapFile);

print STDERR "writing matrix ...\n" if($verbose);
my $newHeaderMatrixFile=$output.".newHeaders.matrix.gz";
writeMatrix($matrix,$new_inc2header,$newHeaderMatrixFile,$missingValue,$commentLine);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);
