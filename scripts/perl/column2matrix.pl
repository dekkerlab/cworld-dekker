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

sub check_options {
    my $opts = shift;
    
    my ($inputList,$verbose,$output,$optionalYHeaderFile,$optionalXHeaderFile,$missingValue);
 
    if( exists($opts->{ inputList }) ) {
        $inputList = $opts->{ inputList };
    } else {
        print STDERR "\nERROR: Option inputList|i is required.\n";
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
    
    if( exists($opts->{ optionalYHeaderFile }) ) {
        $optionalYHeaderFile = $opts->{ optionalYHeaderFile };
    } else {
        $optionalYHeaderFile="";
    }
    
    if( exists($opts->{ optionalXHeaderFile }) ) {
        $optionalXHeaderFile = $opts->{ optionalXHeaderFile };
    } else {
        $optionalXHeaderFile="";
    }
    
    if( exists($opts->{ missingValue }) ) {
        $missingValue = $opts->{ missingValue };
    } else {
        $missingValue="NA";
    }
    
    return($inputList,$verbose,$output,$optionalYHeaderFile,$optionalXHeaderFile,$missingValue);
    
}


sub intro() {
    print STDERR "\n";
    
    print STDERR "Tool:\t\tcolumn2matrix.pl\n";
    print STDERR "Version:\t".$cworld::dekker::VERSION."\n";
    print STDERR "Summary:\tturn list (3 tab) file into matrix\n";
    
    print STDERR "\n";
}

sub help() {
    intro();
    
    print STDERR "Usage: perl column2matrix.pl [OPTIONS] -i <inputMatrix>\n";
    
    print STDERR "\n";
    
    print STDERR "Required:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-i", "[]", "input list file");
    
    print STDERR "\n";
    
    print STDERR "Options:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-v", "[]", "FLAG, verbose mode");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-o", "[]", "prefix for output file(s)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--oxh", "[]", "optional x-axis header file");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--oyh", "[]", "optional y-axis header file");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--mv", "[]", "optional missing-value [0,NA]");

    print STDERR "\n";
    
    print STDERR "Notes:";
    print STDERR "
    This script converts a 3 column txt (tsv) file to a matrix.
    Use -oxh and -oyh to set the x-axis / y-axis headers.
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

sub readOptionalHeaderFile($;$) {
    # required
    my $optionalHeaderFile=shift;
    # optional
    my $verbose=0;
    $verbose=shift if @_;
    
    print STDERR "processing optional header file ($optionalHeaderFile)\n" if($verbose);
    croak "optionalHeaderFile [$optionalHeaderFile] does not exist" if(!(-e $optionalHeaderFile));
    
    my @optionalHeaders=();
    
    open(IN,inputWrapper($optionalHeaderFile)) or croak "Could not open file [$optionalHeaderFile] - $!";
    while(my $line = <IN>) {
        chomp($line);
        next if(($line eq "") or ($line =~ m/^#/));
        
        push(@optionalHeaders,$line);
    }
    close(IN);
    
    my $nOptionalHeaders=@optionalHeaders;
    
    print STDERR "\tfound $nOptionalHeaders\n" if($verbose);
    print STDERR "\n" if($verbose);
    
    return(@optionalHeaders);
}

my %options;
my $results = GetOptions( \%options,'inputList|i=s','verbose|v','output|o=s','optionalYHeaderFile|oyh=s','optionalXHeaderFile|oxh=s','missingValue|mv=s') or croak help();

my ($inputList,$verbose,$output,$optionalYHeaderFile,$optionalXHeaderFile,$missingValue)=check_options( \%options );

intro() if($verbose);

#get the absolute path of the script being used
my $cwd = getcwd();
my $fullScriptPath=abs_path($0);
my @fullScriptPathArr=split(/\//,$fullScriptPath);
@fullScriptPathArr=@fullScriptPathArr[0..@fullScriptPathArr-3];
my $scriptPath=join("/",@fullScriptPathArr);

croak "inputList [$inputList] does not exist" if(!(-e $inputList));

$output=getFileName($inputList) if($output eq "");

my (@yHeadersArr,@xHeadersArr);
@yHeadersArr=readOptionalHeaderFile($optionalYHeaderFile) if($optionalYHeaderFile ne "");
@xHeadersArr=readOptionalHeaderFile($optionalXHeaderFile) if($optionalYHeaderFile ne "");

my %yHeaders=();
my %xHeaders=();
my %encountered=();

my %matrix=();

my $line="";
my $lineNum=1;
my $nHits=0;

print "reading inputList ...\n" if($verbose);

open(IN,inputWrapper($inputList)) or croak "Could not open file [$inputList] - $!";
while($line = <IN>) {
    chomp($line);
    next if(($line eq "") or ($line =~ m/^#/));
    
    my @tmp=split(/\t/,$line);
    croak "bad input file format! [$line]" if(@tmp > 3);
    
    my $yHeader=$tmp[0];
    my $xHeader=$tmp[1];
    my $score=1;
    $score=$tmp[2] if(defined($tmp[2]));
    $score = "NA" if($score !~ (/^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/));
    
    my $key=$yHeader."___".$xHeader;
    croak "found a duplicate! [lineNum=".$lineNum." | ".$key."]" if(exists($encountered{$key}));
    $encountered{$key}=1;
    
    $nHits++ if(($score ne "NA") and ($score != 0));
    
    # if user did not supply an optional header file - then pull headers from the input file
    if($optionalYHeaderFile eq "") {
        if(!exists($yHeaders{$yHeader})) {
            push(@yHeadersArr,$yHeader);
            $yHeaders{$yHeader}=1;
        }
    }
    
    # if user did not supply an optional header file - then pull headers from the input file
    if($optionalXHeaderFile eq "") {
        if(!exists($xHeaders{$xHeader})) {
            push(@xHeadersArr,$xHeader);
            $xHeaders{$xHeader}=1;
        }
    }
    
    $matrix{$yHeader}{$xHeader}=$score;
    
    $lineNum++;
    
}
close(IN);

print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

print STDERR "found $nHits binary [!=0] interactions...\n" if($verbose);

my $nYHeaders=@yHeadersArr;
my $nXHeaders=@xHeadersArr;

print STDERR "found $nYHeaders yHeaders\n" if($verbose);
print STDERR "found $nXHeaders xHeaders\n" if($verbose);

print STDERR "\n" if($verbose);

print STDERR "writing output matrices ...\n" if($verbose);

open(BINARY_MATRIX,outputWrapper($output.".binary.matrix.gz")) or croak "Could not open file [$output] - $!";
open(SCORE_MATRIX,outputWrapper($output.".score.matrix.gz")) or croak "Could not open file [$output] - $!";

my $time=getDate();

print BINARY_MATRIX "# UMASS 5C\n";
print BINARY_MATRIX "# http://3DG.umassmed.edu\n";
print BINARY_MATRIX "# $time\n";

print SCORE_MATRIX "# UMASS 5C\n";
print SCORE_MATRIX "# http://3DG.umassmed.edu\n";
print SCORE_MATRIX "# $time\n";


print BINARY_MATRIX "\t";
print SCORE_MATRIX "\t";
for(my $i=0;$i<$nXHeaders;$i++) {
    print BINARY_MATRIX $xHeadersArr[$i] . "\t";
    print SCORE_MATRIX $xHeadersArr[$i] . "\t";
}
print BINARY_MATRIX "\n";
print SCORE_MATRIX "\n";

for(my $y=0;$y<$nYHeaders;$y++) { 
    my $yHeader=$yHeadersArr[$y];
    print BINARY_MATRIX $yHeadersArr[$y] . "\t";
    print SCORE_MATRIX $yHeadersArr[$y] . "\t";
    for(my $x=0;$x<$nXHeaders;$x++) { 
        my $xHeader=$xHeadersArr[$x];
        
        my $score=$missingValue;
        my $binary="NA";
        $score=$matrix{$yHeader}{$xHeader} if(exists($matrix{$yHeader}{$xHeader}));
        $binary=1 if(exists($matrix{$yHeader}{$xHeader}));
        
        print BINARY_MATRIX $binary;
        print SCORE_MATRIX $score;
        if($x !=($nXHeaders-1)) { print BINARY_MATRIX "\t"; print SCORE_MATRIX "\t"; }
    }
    if($y != ($nYHeaders-1)) { print BINARY_MATRIX "\n"; print SCORE_MATRIX "\n"; }
}
close(BINARY_MATRIX);
close(SCORE_MATRIX);

print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);