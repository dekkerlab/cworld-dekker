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

    my ($inputMatrix,$verbose,$output,$minDistance,$maxDistance,$excludeCis,$excludeTrans,$elementBedFiles,$y_elementBedFiles,$x_elementBedFiles,$elementName,$zoomCoordinates,$y_zoomCoordinates,$x_zoomCoordinates,$elementZoneSize);

        
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
    
    if( exists($opts->{ minDistance }) ) {
        $minDistance = $opts->{ minDistance };
    } else {
        $minDistance=undef;
    }
    
    if( exists($opts->{ maxDistance }) ) {
        $maxDistance = $opts->{ maxDistance };
    } else {
        $maxDistance = undef;
    }
    
    if( exists($opts->{ excludeCis }) ) {
        $excludeCis=1;
    } else {
        $excludeCis=0;
    }
    
    if( exists($opts->{ excludeTrans }) ) {
        $excludeTrans=1;
    } else {
        $excludeTrans=0;
    }
    
    if( exists($opts->{ elementBedFiles }) ) {
        $elementBedFiles = $opts->{ elementBedFiles };
    } else {
        $elementBedFiles = [];
    }
    
    if( exists($opts->{ y_elementBedFiles }) ) {
        $y_elementBedFiles = $opts->{ y_elementBedFiles };
    } else {
        $y_elementBedFiles = [];
    }
    
    if( exists($opts->{ x_elementBedFiles }) ) {
        $x_elementBedFiles = $opts->{ x_elementBedFiles };
    } else {
        $x_elementBedFiles = [];
    }
    
    if( exists($opts->{ elementName }) ) {
        $elementName = $opts->{ elementName };
    } else {
        $elementName=getSmallUniqueString();
    }
    
    if( exists($opts->{ zoomCoordinates }) ) {
        $zoomCoordinates = $opts->{ zoomCoordinates };
    } else {
        $zoomCoordinates = [];
    }
    
    if( exists($opts->{ y_zoomCoordinates }) ) {
        $y_zoomCoordinates = $opts->{ y_zoomCoordinates };
    } else {
        $y_zoomCoordinates = [];
    }
    
    if( exists($opts->{ x_zoomCoordinates }) ) {
        $x_zoomCoordinates = $opts->{ x_zoomCoordinates };
    } else {
        $x_zoomCoordinates = [];
    }
    
    if( exists($opts->{ elementZoneSize }) ) {
        $elementZoneSize = $opts->{ elementZoneSize };
    } else {
        $elementZoneSize=0;
    }
    
    return($inputMatrix,$verbose,$output,$minDistance,$maxDistance,$excludeCis,$excludeTrans,$elementBedFiles,$y_elementBedFiles,$x_elementBedFiles,$elementName,$zoomCoordinates,$y_zoomCoordinates,$x_zoomCoordinates,$elementZoneSize);
}


sub intro() {
    print STDERR "\n";
    
    print STDERR "Tool:\t\tsubsetMatrix.pl\n";
    print STDERR "Version:\t".$cworld::dekker::VERSION."\n";
    print STDERR "Summary:\tsubset matrix by distance, or by BED file (bin overlap)\n";
    
    print STDERR "\n";
}

sub help() {
    intro();
    
    print STDERR "Usage: perl subsetMatrix.pl [OPTIONS] -i <inputMatrix>\n";
    
    print STDERR "\n";
    
    print STDERR "Required:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-i", "[]", "input matrix file");
    
    print STDERR "\n";
    
    print STDERR "Options:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-v", "[]", "FLAG, verbose mode");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-o", "[]", "prefix for output file(s)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--minDist", "[]", "minimum allowed interaction distance, exclude < N distance (in BP)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--maxDist", "[]", "maximum allowed interaction distance, exclude > N distance (in BP)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ec", "[]", "FLAG, exclude CIS data");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--et", "[]", "FLAG, exclude TRANS data");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ebf@", "[]", "x/y axis element bed file");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--yebf@", "[]", "y axis element bed file");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--xebf@", "[]", "x axis element bed file");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--z@", "[]", "x/y axis zoom coordinate [UCSC]");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--yz@", "[]", "y axis zoom coordinate [UCSC]");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--xz@", "[]", "x axis zoom coordinate [UCSC]");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--en", "[]", "elementName, descriptor for output files");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-zs", "[]", "increase element size by N bp, (increase overlap)");
    
    print STDERR "\n";
    
    print STDERR "Notes:";
    print STDERR "
    This script aggregrates cData surrounding a list of 'elements'.
    Input Matrix can be TXT or gzipped TXT.
    AnchorListFile must be BED5+ format.
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

sub processZoom($$$) {
    my $zoomCoordinates=shift;
    my $subset_inc2header=shift;
    my $subset_header2inc=shift;
    
    my $numSubsetHeaders=keys(%{$subset_inc2header});
    
    my $zoom_header2inc={};
    my $zoom_inc2header={};
    
    my $num_zoomHeaders=0;
    
    if(@$zoomCoordinates > 0) {
        for(my $z=0;$z<@$zoomCoordinates;$z++) {
            my $tmp_zoomCoordinate=$zoomCoordinates->[$z];
            my $tmp_zoomObject=splitCoordinate($tmp_zoomCoordinate);
            my $zoomChromosome=$tmp_zoomObject->{ chromosome };
            my $zoomStart=$tmp_zoomObject->{ start };
            my $zoomEnd=$tmp_zoomObject->{ end };
            my $flag=$tmp_zoomObject->{ flag };
            my $name=$tmp_zoomObject->{ name };
            
            for(my $h=0;$h<$numSubsetHeaders;$h++) {
                my $header=$subset_inc2header->{$h};
                my $headerObject=getHeaderObject($header);
                my $headerStart=$headerObject->{ start };
                my $headerEnd=$headerObject->{ end };
                my $headerChromosome=$headerObject->{ chromosome };
                my $stripped_xheaderChromosome=stripChromosomeGroup($headerChromosome);
                                
                next if(($headerChromosome ne $zoomChromosome) and ($stripped_xheaderChromosome ne $zoomChromosome));
                next if($headerEnd < $zoomStart);
                next if($headerStart > $zoomEnd);
                
                if(!exists($zoom_header2inc->{$header})) {
                    $zoom_header2inc->{$header}=$num_zoomHeaders;
                    $zoom_inc2header->{$num_zoomHeaders}=$header;
                    $num_zoomHeaders++;
                }
            }
        }
    }
    
    return($zoom_inc2header,$zoom_header2inc)
    
}

sub processElements($$$$$$) {
    my $elementBedFiles=shift;
    my $elementName=shift;
    my $elementZoneSize=shift;
    my $headerBedFile=shift;
    my $subset_inc2header=shift;
    my $subset_header2inc=shift;
    
    for(my $i=0;$i<@{$elementBedFiles};$i++) {
        my $elementBedFile = $elementBedFiles->[$i];
        print STDERR "validating $elementBedFile ...\n";
        validateBED($elementBedFile);
    }

    print STDERR "\n";

    my $elementBedFile=$elementBedFiles->[0] if(@{$elementBedFiles});
    $elementBedFile=combineBedFiles($elementBedFiles,$elementName) if(@{$elementBedFiles} > 1);

    print STDERR "intersecting BED files ...\n";
    my $bedOverlapFile=intersectBED($headerBedFile,$elementBedFile,$elementZoneSize);
    print STDERR "\t$bedOverlapFile\n";
    system("rm '".$elementBedFile."'") if(@{$elementBedFiles} > 1);

    print STDERR "\n";

    print STDERR "loading BED file ...\n";
    my ($elements)=loadBED($bedOverlapFile);
    print STDERR "\tfound ".@{$elements}." elements\n";
    system("rm '".$bedOverlapFile."'");

    croak "found no overlapping headers!" if(@{$elements} == 0);
    
    my $numElements=@{$elements};
    my $numHeaders=0;
    for(my $e=0;$e<$numElements;$e++) {
        my $elementHeader=$elements->[$e]->{ name };
        next if(exists($subset_header2inc->{$elementHeader}));
        my $elementObject=getHeaderObject($elementHeader,1);
        $subset_inc2header->{$numHeaders}=$elementHeader;
        $subset_header2inc->{$elementHeader}=$numHeaders;
        $numHeaders++;
    }
    
    return($subset_inc2header,$subset_header2inc);

}

my %options;
my $results = GetOptions( \%options,'inputMatrix|i=s','verbose|v','output|o=s','minDistance|minDist=i','maxDistance|maxDist=i','excludeCis|ec','excludeTrans|et','elementBedFiles|ebf=s@','y_elementBedFiles|yebf=s@','x_elementBedFiles|xebf=s@','elementName|en=s','zoomCoordinates|z=s@','y_zoomCoordinates|yz=s@','x_zoomCoordinates|xz=s@','elementZoneSize|ezs=s') or croak help();

my ($inputMatrix,$verbose,$output,$minDistance,$maxDistance,$excludeCis,$excludeTrans,$elementBedFiles,$y_elementBedFiles,$x_elementBedFiles,$elementName,$zoomCoordinates,$y_zoomCoordinates,$x_zoomCoordinates,$elementZoneSize)=check_options( \%options );

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

my $subset_header2inc=$header2inc;
my $subset_inc2header=$inc2header;

print STDERR "running headers2bed ...\n" if($verbose);
my $headerBedFile=headers2bed($matrixObject);
print STDERR "\t$headerBedFile\n" if($verbose);

print STDERR "\n" if($verbose);

# process element bed files
$y_elementBedFiles = [ @$elementBedFiles, @$y_elementBedFiles ];
$x_elementBedFiles = [ @$elementBedFiles, @$x_elementBedFiles ];

# load bed files
($subset_inc2header->{ y },$subset_header2inc->{ y })=processElements($y_elementBedFiles,$elementName,$elementZoneSize,$headerBedFile,$subset_inc2header->{ y },$subset_header2inc->{ y }) if(@{$y_elementBedFiles} > 0);
($subset_inc2header->{ x },$subset_header2inc->{ x })=processElements($x_elementBedFiles,$elementName,$elementZoneSize,$headerBedFile,$subset_inc2header->{ x },$subset_header2inc->{ x }) if(@{$x_elementBedFiles} > 0);
system("rm '".$headerBedFile."'");

# process zoom coordinates
$y_zoomCoordinates = [ @$zoomCoordinates, @$y_zoomCoordinates ];
$x_zoomCoordinates = [ @$zoomCoordinates, @$x_zoomCoordinates ];

# load zoom coordinates
($subset_inc2header->{ y },$subset_header2inc->{ y })=processZoom($y_zoomCoordinates,$subset_inc2header->{ y },$subset_header2inc->{ y }) if(@{$y_zoomCoordinates} > 0);
($subset_inc2header->{ x },$subset_header2inc->{ x })=processZoom($x_zoomCoordinates,$subset_inc2header->{ x },$subset_header2inc->{ x }) if(@{$x_zoomCoordinates} > 0);

# reset headers if subset
my $num_subsetXHeaders=keys(%{$subset_inc2header->{ x }});
my $num_subsetYHeaders=keys(%{$subset_inc2header->{ y }});
$inc2header->{ x }=$subset_inc2header->{ x } if($num_subsetXHeaders != 0);
$inc2header->{ y }=$subset_inc2header->{ y } if($num_subsetYHeaders != 0);
$header2inc->{ x }=$subset_header2inc->{ x } if($num_subsetXHeaders != 0);
$header2inc->{ y }=$subset_header2inc->{ y } if($num_subsetYHeaders != 0);

# update matrix object
$matrixObject->{ inc2header }=$inc2header;
$matrixObject->{ header2inc }=$header2inc;
$matrixObject=updateMatrixObject($matrixObject);
$inc2header=$matrixObject->{ inc2header };
$header2inc=$matrixObject->{ header2inc };
$numYHeaders=$matrixObject->{ numYHeaders };
$numXHeaders=$matrixObject->{ numXHeaders };

print STDERR "\n" if($verbose);

$output .= "__".$elementName;

#read Matrix
my $matrix={};
($matrix)=getData($inputMatrix,$matrixObject,$verbose,$minDistance,$maxDistance,$excludeCis,$excludeTrans);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

my $subsetMatrixFile=$output.".subset.matrix.gz";
$subsetMatrixFile=$output.".subset".$maxDistance.".matrix.gz" if(defined($maxDistance));
print STDERR "Writing matrix to file ($subsetMatrixFile)...\n" if($verbose);
writeMatrix($matrix,$inc2header,$subsetMatrixFile,"NA");
print STDERR "\tcomplete\n\n" if($verbose);

