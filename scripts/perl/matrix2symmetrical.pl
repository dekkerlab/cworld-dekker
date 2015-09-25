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

    my ($inputMatrix,$verbose,$output);
    
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
    
    $ret->{ inputMatrix }=$inputMatrix;
    $ret->{ verbose }=$verbose;
    $ret->{ output }=$output;
    
    return($ret,$inputMatrix,$verbose,$output);
}


sub intro() {
    print STDERR "\n";
    
    print STDERR "Tool:\t\t".$tool."\n";
    print STDERR "Version:\t".$cworld::dekker::VERSION."\n";
    print STDERR "Summary:\ttransform rectangular matrix into square (symmetrical) matrix\n";
    
    print STDERR "\n";
}

sub help() {
    intro();
    
    print STDERR "Usage: perl matrix2symmetrical.pl [OPTIONS] -i <inputMatrix>\n";
    
    print STDERR "\n";
    
    print STDERR "Required:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-i", "[]", "input matrix file");
    
    print STDERR "\n";
    
    print STDERR "Options:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-v", "[]", "FLAG, verbose mode");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-o", "[]", "prefix for output file(s)");

    print STDERR "\n";
    
    print STDERR "Notes:";
    print STDERR "
    This script turns any matrix into a symmetrical matrix.  If the matrix is already symmetrical, nothing is changed.
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

sub symmetricizeData($$) {
    my $matrixObject=shift;
    my $matrix=shift;
    
    my $inc2header=$matrixObject->{ inc2header };
    my $header2inc=$matrixObject->{ header2inc };
    my $numYHeaders=$matrixObject->{ numYHeaders };
    my $numXHeaders=$matrixObject->{ numXHeaders };    
    my $missingValue=$matrixObject->{ missingValue };
    
    for(my $y=0;$y<$numYHeaders;$y++) {
        my $yHeader=$inc2header->{ y }->{$y};
        for(my $x=0;$x<$numXHeaders;$x++) {
            my $xHeader=$inc2header->{ x }->{$x};    
        
            my $cScore=$matrixObject->{ missingValue };
            $cScore=$matrix->{$y}->{$x} if(defined($matrix->{$y}->{$x}));
            $cScore=$matrix->{$x}->{$y} if(defined($matrix->{$x}->{$y}));
            
            $matrix->{$x}->{$y}=$cScore;
            $matrix->{$y}->{$x}=$cScore;
        }
    }
    
    return($matrix);
}

my %options;
my $results = GetOptions( \%options,'inputMatrix|i=s','verbose|v','output|o=s') or croak help();
my ($ret,$inputMatrix,$verbose,$output)=check_options( \%options );

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
my $numTotalHeaders=$matrixObject->{ numTotalHeaders };
my $missingValue=$matrixObject->{ missingValue };
my $symmetrical=$matrixObject->{ symmetrical };
my $equalHeaderFlag=$matrixObject->{ equalHeaderFlag };
my $inputMatrixName=$matrixObject->{ inputMatrixName };
$output=$matrixObject->{ output };

$inc2header->{ x } = $inc2header->{ xy };
$inc2header->{ y } = $inc2header->{ xy };
$header2inc->{ x } = $header2inc->{ xy };
$header2inc->{ y } = $header2inc->{ xy };
$numYHeaders=$numTotalHeaders;
$numXHeaders=$numTotalHeaders;

$matrixObject->{ inc2header }=$inc2header;
$matrixObject->{ header2inc }=$header2inc;
$matrixObject->{ numYHeaders }=$numYHeaders;
$matrixObject->{ numXHeaders }=$numXHeaders;
$matrixObject->{ missingValue }="NA";

#read Matrix
my $matrix={};
($matrix)=getData($inputMatrix,$matrixObject,$verbose);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

print STDERR "symmetriciz'n data...\n" if($verbose);
$matrix=symmetricizeData($matrixObject,$matrix);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

my $symmetricalMatrixFile = $output.".symmetrical.matrix.gz";
print STDERR "writing matrix ($symmetricalMatrixFile)...\n" if($verbose);
writeMatrix($matrix,$inc2header,$symmetricalMatrixFile,$matrixObject->{ missingValue });
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);