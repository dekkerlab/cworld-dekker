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

    my ($inputMatrixArray,$verbose,$output,$excludeDiagonal);
    
    my $ret={};
    
    if( exists($opts->{ inputMatrixArray }) ) {
        $inputMatrixArray = $opts->{ inputMatrixArray };
    } else {
        print STDERR "\nERROR: Option inputMatrixArray|i is required.\n";
        help()
    }
    
    if( exists($opts->{ verbose }) ) {
        $verbose = $opts->{ verbose };
    } else {
        $verbose = 0;
    }
    
    if( exists($opts->{ output }) ) {
        $output = $opts->{ output };
    } else {
        $output = "default";
    }
    
    if( exists($opts->{ excludeDiagonal }) ) {
        $excludeDiagonal = 1;
    } else {
        $excludeDiagonal = 0;
    }
    
    $ret->{ inputMatrixArray }=$inputMatrixArray;
    $ret->{ verbose }=$verbose;
    $ret->{ output }=$output;
    $ret->{ excludeDiagonal }=$excludeDiagonal;
    
    return($ret,$inputMatrixArray,$verbose,$output,$excludeDiagonal);
}

sub intro() {
    print STDERR "\n";
    
    print STDERR "Tool:\t\t".$tool."\n";
    print STDERR "Version:\t".$cworld::dekker::VERSION."\n";
    print STDERR "Summary:\tnormalizes matrix sum - scales to 10^6\n";
    
    print STDERR "\n";
}

sub help() {
    intro();
    
    print STDERR "Usage: perl normalizeMatrix.pl [OPTIONS] -i <inputMatrix>\n";
    
    print STDERR "\n";
    
    print STDERR "Required:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-i", "[]", "input matrix file, multiple files allowed");
    
    print STDERR "\n";
    
    print STDERR "Options:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-v", "[]", "FLAG, verbose mode");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-o", "[]", "prefix for output file(s)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--id", "[]", "FLAG, ignore diagonal bin during normalization");
    
    print STDERR "\n";
    
    print STDERR "Notes:";
    print STDERR "
    This script normalizes the sum of a matrix.
    If the matrix is symmetrical, then only half of the matrix is summed (diagonal + x>y data)
    Each interaction is first translated into a percentage of the current sum.
    Each value is then multipled by 1,000,000 so that the new sum is 10^6
    Matrix should have row/col headers in the standard my5C format.
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

sub getNAs($$$) {
    my $matrixObject=shift;
    my $matrix=shift;
    my $na_matrix=shift;
    
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
            
            $na_matrix->{$y}->{$x}=1 if($cScore eq "NA");
            
        }
    }
    
    return($na_matrix);
}

sub applyNAMask($$$) {
    my $matrixObject=shift;
    my $matrix=shift;
    my $na_matrix=shift;
    
    my $inc2header=$matrixObject->{ inc2header };
    my $header2inc=$matrixObject->{ header2inc };
    my $numYHeaders=$matrixObject->{ numYHeaders };
    my $numXHeaders=$matrixObject->{ numXHeaders };
    my $missingValue=$matrixObject->{ missingValue };
    
    for(my $y=0;$y<$numYHeaders;$y++) {
        my $yHeader=$inc2header->{ y }->{$y};        
        for(my $x=0;$x<$numXHeaders;$x++) {
            my $xHeader=$inc2header->{ x }->{$x};    
        
            $matrix->{$y}->{$x}="NA" if(exists($na_matrix->{$y}{$x}));
            delete($matrix->{$y}->{$x}) if((exists($na_matrix->{$y}{$x})) and ($missingValue eq "NA"));
        }
    }
    
    return($matrix);
}


my %options;
my $results = GetOptions( \%options,'inputMatrixArray|i=s@','verbose|v','output|o=s','excludeDiagonal|ed') or croak help();
my ($ret,$inputMatrixArray,$verbose,$output,$excludeDiagonal)=check_options( \%options );

intro() if($verbose);

#get the absolute path of the script being used
my $cwd = getcwd();
my $fullScriptPath=abs_path($0);
my @fullScriptPathArr=split(/\//,$fullScriptPath);
@fullScriptPathArr=@fullScriptPathArr[0..@fullScriptPathArr-3];
my $scriptPath=join("/",@fullScriptPathArr);
my $commentLine=getScriptOpts($ret,$tool);

my $na_matrix={};

if(@{$inputMatrixArray} > 1) {
    
    # build the NA mask 
    print STDERR "building NA mask\n" if($verbose);
    for(my $i=0;$i<@{$inputMatrixArray};$i++) {
        my $inputMatrix = $inputMatrixArray->[$i];
         
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

        my $matrix={};
        ($matrix)=getData($inputMatrix,$matrixObject);
        $na_matrix=getNAs($matrixObject,$matrix,$na_matrix);
        
    }

    print STDERR "\n";
}

#since all matrices are same size, use info from first
my $matrixObject=getMatrixObject($inputMatrixArray->[0],$output,$verbose);
my $inc2header=$matrixObject->{ inc2header };
my $header2inc=$matrixObject->{ header2inc };
my $numYHeaders=$matrixObject->{ numYHeaders };
my $numXHeaders=$matrixObject->{ numXHeaders };

if(@{$inputMatrixArray} > 1) {
    my $na_matrixFile=$output.".nan.matrix.gz";
    writeMatrix($na_matrix,$inc2header,$na_matrixFile,"NA",$commentLine);
    my $nNAs=getMatrixSum($matrixObject,$na_matrix);
    print STDERR "NA mask\n" if($verbose);
    print STDERR "\t$nNAs NAs\n" if($verbose);

    print STDERR "\n" if($verbose);
}

if(@{$inputMatrixArray} > 1) {
    # now apply the NA mask to all matrices
    print STDERR "applying NA mask...\n" if($verbose);
    for(my $i=0;$i<@{$inputMatrixArray};$i++) {
        my $inputMatrix = $inputMatrixArray->[$i];
        
        croak "inputMatrix [$inputMatrix] does not exist" if(!(-e $inputMatrix));

        # get matrix information
        my $matrixObject=getMatrixObject($inputMatrix,"",0);
        my $inc2header=$matrixObject->{ inc2header };
        my $header2inc=$matrixObject->{ header2inc };
        my $numYHeaders=$matrixObject->{ numYHeaders };
        my $numXHeaders=$matrixObject->{ numXHeaders };
        my $missingValue=$matrixObject->{ missingValue };
        my $symmetrical=$matrixObject->{ symmetrical };
        my $equalHeaderFlag=$matrixObject->{ equalHeaderFlag };
        my $inputMatrixName=$matrixObject->{ inputMatrixName };
        $output=$matrixObject->{ output };

        my ($matrix)=getData($inputMatrix,$matrixObject);
        $matrix=applyNAMask($matrixObject,$matrix,$na_matrix);

        my $maskedFile=$output.".masked.matrix.gz";
        writeMatrix($matrix,$inc2header,$maskedFile,$missingValue,$commentLine);
        # replace input file with masked file
        $inputMatrixArray->[$i]=$maskedFile;
        
        undef $matrix;
        
    }

    print STDERR "\n" if($verbose);

    # validate all files are the same
    print STDERR "validating identical matrices...\n" if($verbose);
    for(my $i1=0;$i1<@{$inputMatrixArray};$i1++) {
        my $inputMatrix_1 = $inputMatrixArray->[$i1];
        for(my $i2=0;$i2<@{$inputMatrixArray};$i2++) {
            my $inputMatrix_2 = $inputMatrixArray->[$i2];
            
            next if($inputMatrix_1 eq $inputMatrix_2);
            
            validateIdenticalMatrixStructure($inputMatrix_1,$inputMatrix_2);
        }
    }
    print STDERR "\tdone\n" if($verbose);

    print STDERR "\n" if($verbose);
}

print STDERR "normalizing...\n" if($verbose);
my @normalizedMatrixArray=();
for(my $i=0;$i<@{$inputMatrixArray};$i++) {
    my $inputMatrix = $inputMatrixArray->[$i];
    
    croak "inputMatrix [$inputMatrix] does not exist" if(!(-e $inputMatrix));

    # get matrix information
    my $matrixObject=getMatrixObject($inputMatrix,"",0);
    my $inc2header=$matrixObject->{ inc2header };
    my $header2inc=$matrixObject->{ header2inc };
    my $numYHeaders=$matrixObject->{ numYHeaders };
    my $numXHeaders=$matrixObject->{ numXHeaders };
    my $missingValue=$matrixObject->{ missingValue };
    my $symmetrical=$matrixObject->{ symmetrical };
    my $equalHeaderFlag=$matrixObject->{ equalHeaderFlag };
    my $inputMatrixName=$matrixObject->{ inputMatrixName };
    $output=$matrixObject->{ output };

    #read Matrix
    my ($matrix)=getData($inputMatrix,$matrixObject,$verbose);
    
    #normalize the matrix
    $matrix=normalizeMatrix($matrixObject,$matrix,$excludeDiagonal);
    
    #write the normalized matrix
    my $normalizedMatrixFile=$output.".normalized.matrix.gz";
    writeMatrix($matrix,$inc2header,$normalizedMatrixFile,$missingValue,$commentLine);
    
    undef $matrix;
    
    push(@normalizedMatrixArray,$normalizedMatrixFile);

}

print STDERR "\n" if($verbose);

for(my $i=0;$i<@normalizedMatrixArray;$i++) {
    my $inputMatrix = $normalizedMatrixArray[$i];
    
    croak "inputMatrix [$inputMatrix] does not exist" if(!(-e $inputMatrix));

    # get matrix information
    my $matrixObject=getMatrixObject($inputMatrix,"",0);
    my $inc2header=$matrixObject->{ inc2header };
    my $header2inc=$matrixObject->{ header2inc };
    my $numYHeaders=$matrixObject->{ numYHeaders };
    my $numXHeaders=$matrixObject->{ numXHeaders };
    my $missingValue=$matrixObject->{ missingValue };
    my $symmetrical=$matrixObject->{ symmetrical };
    my $equalHeaderFlag=$matrixObject->{ equalHeaderFlag };
    my $inputMatrixName=$matrixObject->{ inputMatrixName };
    $output=$matrixObject->{ output };
    
    print STDERR "$output\n" if($verbose);
    
    if($numYHeaders < 2000) {
        my ($matrix)=getData($inputMatrix,$matrixObject);        
        my ($matrixSum)=getMatrixSum($matrixObject,$matrix,$excludeDiagonal);
    
        print STDERR "\tmatrixSum\t$matrixSum\n" if($verbose);
    
        my ($totalReads,$cisPercent,$transPercent,$averageTrans) = getMatrixAttributes($inputMatrix);
        print STDERR "\tcis\t".$cisPercent."%\n" if($verbose);
        print STDERR "\ttrans\t".$transPercent."%\n" if($verbose);
        print STDERR "\taverageTransSignal\t$averageTrans\n" if($verbose);
        print STDERR "\n" if($verbose);
    }
}

print STDERR "\n" if($verbose);