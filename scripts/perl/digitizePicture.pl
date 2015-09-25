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

use GD::Simple;

use cworld::dekker;

my $tool=(split(/\//,abs_path($0)))[-1];

sub check_options {
    my $opts = shift;
    
    my ($inputPNG,$verbose,$output,$colorMode);
    
    my $ret={};
    
    if( exists($opts->{ inputPNG }) ) {
        $inputPNG = $opts->{ inputPNG };
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
    
    if( exists($opts->{ colorMode }) ) {
        $colorMode = $opts->{ colorMode };
    } else {
        $colorMode="mean";
    }
    
    $ret->{ inputPNG }=$inputPNG;
    $ret->{ verbose }=$verbose;
    $ret->{ output }=$output;
    $ret->{ colorMode }=$colorMode;
    
    return($ret,$inputPNG,$verbose,$output,$colorMode);
}

sub PNG2matrix($$) {
    my $png=shift;
    my $colorMode=shift;
    
    my $image = GD::Image->newFromPng($png, 1); 
    my $width = $image->width;
    my $height = $image->height;

    my %matrix=();
    my ($x,$y);
    for($y=0;$y<$height;$y++) {
        for($x=0;$x<$width;$x++) {
            my $index = $image->getPixel($x, $y);
            my ($red,$green,$blue) = $image->rgb($index);
            
            $red = 255 - $red;
            $green = 255 - $green;
            $blue = 255 - $blue;
            
            my $color=int(($red+$green+$blue)/3);
            
            my @tmp=($red,$green,$blue);
            @tmp = sort {$a <=> $b} @tmp;
            
            $color=$tmp[1] if($colorMode eq "median");
            $color=int(($red+$green+$blue)/3) if($colorMode eq "mean");
            $color=$red if($colorMode eq "red");
            $color=$green if($colorMode eq "green");
            $color=$blue if($colorMode eq "blue");
            $color=max($red,$green,$blue) if($colorMode eq "max");
            $color=min($red,$green,$blue) if($colorMode eq "min");
            $color=min($red,$green,$blue) if($colorMode eq "median");
            
            $matrix{$y}{$x}=$color;
        }
    }
    
    return($height,$width,\%matrix);
}

sub matrix2TXT($$$$) {
    my $matrix=shift;
    my $height=shift;
    my $width=shift;
    my $outFile=shift;
    
    open(OUT,outputWrapper($outFile)) or croak "Could not open file [$outFile] - $!";
    
    for(my $x=0;$x<$width;$x++) {
        my $xHeader="x".$x."|NA|NA:".$x."-".($x+1);
        
        print OUT "\t".$xHeader;
    }
    
    print OUT "\n";

    for(my $y=0;$y<$height;$y++) {
        
        my $yHeader="y".$y."|NA|NA:".$y."-".($y+1);
        print OUT $yHeader;
        
        for(my $x=0;$x<$width;$x++) {
            my $xHeader="x".$x."|NA|NA:".$x."-".($x+1);
            
            my $value = "NA";
            $value = $value=$matrix->{$y}->{$x} if(exists($matrix->{$y}->{$x}));
            
            print OUT "\t$value";
        }
        
        print OUT "\n" if($y != ($height-1));
    }
    
    close(OUT);
}

sub intro() {
    print STDERR "\n";
    
    print STDERR "Tool:\t\t".$tool."\n";
    print STDERR "Version:\t".$cworld::dekker::VERSION."\n";
    print STDERR "Summary:\tdigitize picture into my5C matrix format\n";
    
    print STDERR "\n";
}

sub help() {
    intro();
    
    print STDERR "Usage: perl digitizePicture.pl [OPTIONS] -i <inputPNG>\n";
    
    print STDERR "\n";
    
    print STDERR "Required:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-i", "[]", "input matrix file");
    
    print STDERR "\n";
    
    print STDERR "Options:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-v", "[]", "FLAG, verbose mode");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-o", "[]", "prefix for output file(s)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--cm", "[mean]", "color mode (red,blue,green,min,max,mean,median)");
    print STDERR "\n";
    
    print STDERR "Notes:";
    print STDERR "
    This script turns any PNG into a my5C formatted matrix.\n";
    
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
my $results = GetOptions( \%options,'inputPNG|i=s','verbose|v','output|o=s','colorMode|cm=s') or croak help();
my ($ret,$inputPNG,$verbose,$output,$colorMode)=check_options( \%options );

intro() if($verbose);

#get the absolute path of the script being used
my $cwd = getcwd();
my $fullScriptPath=abs_path($0);
my @fullScriptPathArr=split(/\//,$fullScriptPath);
@fullScriptPathArr=@fullScriptPathArr[0..@fullScriptPathArr-3];
my $scriptPath=join("/",@fullScriptPathArr);

croak "inputPNG [$inputPNG] does not exist" if(!(-e $inputPNG));
croak "must use PNG image" if($inputPNG !~ /.png$/);

$output=getFileName($inputPNG) if($output eq "");

#read the PNG image in
print STDERR "reading in input image ($inputPNG)...\n" if($verbose);
my ($height,$width,$matrix)=PNG2matrix($inputPNG,$colorMode);
print STDERR "\theight\t$height\n" if($verbose);
print STDERR "\twidth\t$width\n" if($verbose);

print STDERR "\n" if($verbose);

print STDERR "writing matrix ...\n" if($verbose);
matrix2TXT($matrix,$height,$width,$output.".".$colorMode.".digitized.matrix");
print STDERR "\tdone.\n" if($verbose);

print STDERR "\n" if($verbose);