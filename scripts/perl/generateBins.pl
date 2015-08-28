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

    my ($assembly,$verbose,$output,$regionCoordinates,$binSize,$binStep);

    if( exists($opts->{ assembly }) ) {
        $assembly = $opts->{ assembly };
    } else {
        print STDERR "\nERROR: Option assembly|a is required.\n";
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
        $output="myHeaders";
    }
    
    if( exists($opts->{ regionCoordinates }) ) {
        $regionCoordinates = $opts->{ regionCoordinates };
    } else {
        print STDERR "\nERROR: Option regionCoordinates|r is required.\n";
        help();
    }
    
    if( exists($opts->{ binSize }) ) {
        $binSize = $opts->{ binSize };
    } else {
        print STDERR "\nERROR: Option binSize|bsize is required.\n";
        help();
    }
    
    if( exists($opts->{ binStep }) ) {
        $binStep = $opts->{ binStep };
    } else {
        print STDERR "\nERROR: Option binStep|bstep is required.\n";
        help();
    }

    return($assembly,$verbose,$output,$regionCoordinates,$binSize,$binStep);
}


sub intro() {
    print STDERR "\n";
    
    print STDERR "Tool:\t\tgenerateBins.pl\n";
    print STDERR "Version:\t".$cworld::dekker::VERSION."\n";
    print STDERR "Summary:\tcreate my5C formatted headers\n";
    
    print STDERR "\n";
}

sub help() {
    intro();
    
    print STDERR "Usage: perl generateBins.pl [OPTIONS] -a <assembly> -r <regionCoordinates -bsize <binSize> -bstep <binStep>\n";
    
    print STDERR "\n";
    
    print STDERR "Required:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "--a", "[1]", "genome assembly");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--r", "[]", "region coordinates, e.g. chr1:1-50000000 (UCSC)");    
    printf STDERR ("\t%-10s %-10s %-10s\n", "--bsize", "[]", "bin size in bp");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--bstep", "[]", "bin step in factor of overlap, 1=non overlapping, 2-inf=overlapping");
    
    print STDERR "\n";
    
    print STDERR "Options:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-v", "[]", "FLAG, verbose mode");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--n", "[myHeaders]", "optional, prefix for output files");
    
    print STDERR "\n";
    
    print STDERR "Notes:";
    print STDERR "
    This script can create my5C formatted headers.
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

sub generateBins($$$$$) {
    my $assembly=shift;
    my $regionCoordinates=shift;
    my $binSize=shift;
    my $binStep=shift;
    my $output=shift;
    
    my ($coordinateData)=splitCoordinate($regionCoordinates);
    my $chromosome=$coordinateData->{ chromosome };
    my $regionStart=$coordinateData->{ start };
    my $regionEnd=$coordinateData->{ end };


    my $regionLength=($regionEnd-$regionStart);
    my $numBins = ceil($regionLength/($binSize/$binStep));
    
    my $headerBedFile=$output.".headers.bed";
    open(BED,outputWrapper($headerBedFile)) or croak "Could not open file [$headerBedFile] - $!";
    
    print BED "track name='".$output."-headers' description='".$output."-headers ($binSize | $binStep)' maxHeightPixels=128:64:32 visibility=squish autoScale=off useScore=0 color=0,0,0\n";
    
    for(my $i=0;$i<$numBins;$i++) {
        my $binStart = $regionStart+($i*($binSize/$binStep));
        my $binEnd = $binStart + $binSize;
        $binEnd=$regionEnd if($binEnd > $regionEnd);
        my $binIndex=$i;
        my $binName=$binIndex."|".$assembly."|".$chromosome.":".$binStart."-".$binEnd;
        print BED "$chromosome\t$binStart\t$binEnd\t$binName\n";
        
    }
    close(BED);
}

my %options;
my $results = GetOptions( \%options,'assembly|a=s','verbose|v','output|o=s','regionCoordinates|r=s','binSize|bsize=s','binStep|bstep=s') or croak help();

my ($assembly,$verbose,$output,$regionCoordinates,$binSize,$binStep) = check_options( \%options );

intro() if($verbose);

#get the absolute path of the script being used
my $cwd = getcwd();
my $fullScriptPath=abs_path($0);
my @fullScriptPathArr=split(/\//,$fullScriptPath);
@fullScriptPathArr=@fullScriptPathArr[0..@fullScriptPathArr-3];
my $scriptPath=join("/",@fullScriptPathArr);

print STDERR "writing bin positions ... \n" if($verbose);
generateBins($assembly,$regionCoordinates,$binSize,$binStep,$output);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

