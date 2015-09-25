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
    
    my ($inputMatrix,$verbose,$output,$xOrderedHeaderList,$yOrderedHeaderList);
 
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
    
    if( exists($opts->{ xOrderedHeaderList }) ) {
        $xOrderedHeaderList = $opts->{ xOrderedHeaderList };
    } else {
        print STDERR "Option xOrderedHeaderList|xohl is required.\n";
        help();
    }
    
    if( exists($opts->{ yOrderedHeaderList }) ) {
        $yOrderedHeaderList = $opts->{ yOrderedHeaderList };
    } else {
        print STDERR "Option yOrderedHeaderList|yohl is required.\n";
        help();
    }
    
    $ret->{ inputMatrix }=$inputMatrix;
    $ret->{ verbose }=$verbose;
    $ret->{ output }=$output;
    $ret->{ xOrderedHeaderList }=$xOrderedHeaderList;
    $ret->{ yOrderedHeaderList }=$yOrderedHeaderList;
    
    return($ret,$inputMatrix,$verbose,$output,$xOrderedHeaderList,$yOrderedHeaderList);
}

sub intro() {
    print STDERR "\n";
    
    print STDERR "Tool:\t\t".$tool."\n";
    print STDERR "Version:\t".$cworld::dekker::VERSION."\n";
    print STDERR "Summary:\tre-order matrix by list of headers\n";
    
    print STDERR "\n";
}

sub help() {
    intro();
    
    print STDERR "Usage: perl reOrderMatrix.pl [OPTIONS] -i <inputMatrix>\n";
    
    print STDERR "\n";
    
    print STDERR "Required:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-i", "[]", "input matrix file");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--xohl", "[]", "x-axis headers list file");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--yohl", "[]", "y-axis headers list file");
    
    print STDERR "\n";
    
    print STDERR "Options:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-v", "[]", "FLAG, verbose mode");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-o", "[]", "prefix for output file(s)");
    
    print STDERR "\n";
    
    print STDERR "Notes:";
    print STDERR "
    This script can layout primers in 96-well plate format
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

sub getOrderedHeaders($$$) {
    my $headerFile=shift;
    my $header2inc=shift;
    my $axis=shift;
    
    my %orderedHeaders=();
    my %header2OrderedHeader=();
    
    my $nHeaders=0;
    
    open(IN,inputWrapper($headerFile)) or croak "Could not open file [$headerFile] - $!";
    while(my $line = <IN>) {
        chomp($line);
        next if(($line eq "") or ($line =~ m/^#/));
        
        my @tmp=split(/\t/,$line);
        croak "incorrect line format [$line]" if(@tmp > 2);
        
        my $header=$tmp[0];
        my $supplementalHeader=$tmp[0];
        $supplementalHeader=$tmp[1] if(@tmp == 2);
        
        #print STDERR "\t$axis\torderedHeader does not exists in inputMatrix! ($header)\n" if(!(exists($header2inc->{$header})));
        #next if(!(exists($header2inc->{$header})));
        
        $orderedHeaders{$nHeaders}=$header;
        $header2OrderedHeader{$header}=$supplementalHeader;
        
        $nHeaders++;        
    }
    close(IN);
    
    foreach my $header ( keys %{$header2inc}) {
        #print STDERR "\t$axis\tinputMatrix header does not exists in orderedList! ($header)\n" if(!(exists($header2OrderedHeader{$header})));
    }
    
    
    return(\%orderedHeaders,\%header2OrderedHeader);
    
}

my %options;
my $results = GetOptions( \%options,'inputMatrix|i=s','verbose|v','output|o=s','xOrderedHeaderList|xohl=s','yOrderedHeaderList|yohl=s') or croak help();
my ($ret,$inputMatrix,$verbose,$output,$xOrderedHeaderList,$yOrderedHeaderList)=check_options( \%options );

intro() if($verbose);

#get the absolute path of the script being used
my $cwd = getcwd();
my $fullScriptPath=abs_path($0);
my @fullScriptPathArr=split(/\//,$fullScriptPath);
@fullScriptPathArr=@fullScriptPathArr[0..@fullScriptPathArr-3];
my $scriptPath=join("/",@fullScriptPathArr);

croak "inputMatrix [$inputMatrix] does not exist" if(!(-e $inputMatrix));
croak "xOrderedHeaderList [$xOrderedHeaderList] does not exist" if(!(-e $xOrderedHeaderList));
croak "yOrderedHeaderList [$yOrderedHeaderList] does not exist" if(!(-e $yOrderedHeaderList));

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

#read Matrix
my $matrix={};
($matrix)=getData($inputMatrix,$matrixObject,$verbose);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

#read ordered header list
print STDERR "reading ordered headers...\n" if($verbose);
my ($xOrderedHeaders,$xHeader2OrderedHeader)=getOrderedHeaders($xOrderedHeaderList,$header2inc->{ x },'x');
my ($yOrderedHeaders,$yHeader2OrderedHeader)=getOrderedHeaders($yOrderedHeaderList,$header2inc->{ y },'y');
print STDERR "\tdone\n" if($verbose);

my $numOrderedXHeaders=keys %{$xOrderedHeaders};
my $numOrderedYHeaders=keys %{$yOrderedHeaders};

print STDERR "\n" if($verbose);

print STDERR "numYHeaders\t$numYHeaders\t$numOrderedYHeaders\n" if($verbose);
print STDERR "numXHeaders\t$numXHeaders\t$numOrderedXHeaders\n" if($verbose);

print STDERR "\n" if($verbose);

my $reOrderedMatrixFile=$output.".reOrdered.matrix.gz";
open(MATRIX,outputWrapper($reOrderedMatrixFile)) or croak "Could not open file [$reOrderedMatrixFile] - $!";

print STDERR "writing ".$reOrderedMatrixFile." ...\n" if($verbose);

for(my $x=0;$x<$numOrderedXHeaders;$x++) {
    my $xHeader=$xOrderedHeaders->{$x};
    my $xHeaderSupplement=$xHeader2OrderedHeader->{$xHeader};
    print MATRIX "\t$xHeaderSupplement";
}
print MATRIX "\n";

for(my $y=0;$y<$numOrderedYHeaders;$y++) {
    my $yHeader=$yOrderedHeaders->{$y};
    my $yHeaderSupplement=$yHeader2OrderedHeader->{$yHeader};
    
    print MATRIX "$yHeaderSupplement";
    
    for(my $x=0;$x<$numOrderedXHeaders;$x++) {
        my $xHeader=$xOrderedHeaders->{$x};    
        my $xHeaderSupplement=$xHeader2OrderedHeader->{$xHeader};
        
        my $yHeaderIndex=-1;
        $yHeaderIndex=$header2inc->{ y }->{$yHeader} if(exists($header2inc->{ y }->{$yHeader}));
        my $xHeaderIndex=-1;
        $xHeaderIndex=$header2inc->{ x }->{$xHeader} if(exists($header2inc->{ x }->{$xHeader}));
       
        my $cScore=$missingValue;
        $cScore=$matrix->{$yHeaderIndex}->{$xHeaderIndex} if( ($yHeaderIndex != -1) and ($xHeaderIndex != -1) and (exists($matrix->{$yHeaderIndex}->{$xHeaderIndex})) );
        
        $cScore="NA" if( (!exists($header2inc->{ x }->{$xHeader})) or (!exists($header2inc->{ y }->{$yHeader})) );
        $cScore = "NA" if(($cScore eq "") or ($cScore =~ /^NULL$/i) or ($cScore =~ /^NA$/i) or ($cScore =~ /inf$/i) or ($cScore =~ /^nan$/i) or ($cScore eq "NA"));
        
        print MATRIX "\t$cScore";
        
    }
    print MATRIX "\n";
}

close(MATRIX);

print STDERR "\tdone.\n" if($verbose);

print STDERR "\n" if($verbose);