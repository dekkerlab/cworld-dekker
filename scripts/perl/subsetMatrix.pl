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
use Data::Dumper;

use cworld::dekker;

my $tool=(split(/\//,abs_path($0)))[-1];

sub check_options {
    my $opts = shift;

    my ($inputMatrix,$verbose,$output,$minDistance,$maxDistance,$excludeCis,$excludeTrans,$elementBedFiles,$y_elementBedFiles,$x_elementBedFiles,$outputSuffix,$zoomCoordinates,$y_zoomCoordinates,$x_zoomCoordinates,$elementExtension);

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
    
    if( exists($opts->{ outputSuffix }) ) {
        $outputSuffix = $opts->{ outputSuffix };
    } else {
        $outputSuffix = "";
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
    
    if( exists($opts->{ elementExtension }) ) {
        $elementExtension = $opts->{ elementExtension };
    } else {
        $elementExtension=0;
    }
        
    $ret->{ inputMatrix }=$inputMatrix;
    $ret->{ verbose }=$verbose;
    $ret->{ output }=$output;
    $ret->{ minDistance }=$minDistance;
    $ret->{ maxDistance }=$maxDistance;
    $ret->{ excludeCis }=$excludeCis;
    $ret->{ excludeTrans }=$excludeTrans;
    $ret->{ elementBedFiles }=$elementBedFiles;
    $ret->{ y_elementBedFiles }=$y_elementBedFiles;
    $ret->{ x_elementBedFiles }=$x_elementBedFiles;
    $ret->{ outputSuffix }=$outputSuffix;
    $ret->{ zoomCoordinates }=$zoomCoordinates;
    $ret->{ y_zoomCoordinates }=$y_zoomCoordinates;
    $ret->{ x_zoomCoordinates }=$x_zoomCoordinates;
    $ret->{ elementExtension }=$elementExtension;
    
    return($ret,$inputMatrix,$verbose,$output,$minDistance,$maxDistance,$excludeCis,$excludeTrans,$elementBedFiles,$y_elementBedFiles,$x_elementBedFiles,$outputSuffix,$zoomCoordinates,$y_zoomCoordinates,$x_zoomCoordinates,$elementExtension);
}


sub intro() {
    print STDERR "\n";
    
    print STDERR "Tool:\t\t".$tool."\n";
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
    printf STDERR ("\t%-10s %-10s %-10s\n", "--os", "[]", "outputSuffix, suffix for output file");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ee", "[]", "element extensnion, increase element size by N bp, (increase overlap)");
    
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
    my $axis=shift;
    my $elementBedFiles=shift;
    my $outputSuffix=shift;
    my $elementExtension=shift;
    my $headerBedFile=shift;
    my $verbose=shift;
    
    for(my $i=0;$i<@{$elementBedFiles};$i++) {
        my $elementBedFile = $elementBedFiles->[$i];
        print STDERR "validating $elementBedFile ...\n" if($verbose);
        validateBED($elementBedFile);
    }

    print STDERR "\n" if($verbose);

    my $elementBedFile=$elementBedFiles->[0] if(@{$elementBedFiles});
    $elementBedFile=combineBedFiles($elementBedFiles,$outputSuffix) if(@{$elementBedFiles} > 1);

    print STDERR "intersecting BED files ...\n" if($verbose);
    my $bedOverlapFile=intersectBED($headerBedFile,$elementBedFile,$elementExtension);
    print STDERR "\t$bedOverlapFile\n" if($verbose);
    system("rm '".$elementBedFile."'") if(@{$elementBedFiles} > 1);

    print STDERR "\n" if($verbose);

    print STDERR "loading BED file ...\n" if($verbose);
    my ($elements)=loadBED($bedOverlapFile);
    print STDERR "\tfound ".@{$elements}." elements\n" if($verbose);
    system("rm '".$bedOverlapFile."'");

    croak "found no overlapping headers!" if(@{$elements} == 0);
    
    my $element_headers={};
    
    my $numElements=@{$elements};    
    for(my $e=0;$e<$numElements;$e++) {
        my $elementHeader=$elements->[$e]->{ name };
        $element_headers->{$elementHeader}=1;
    }
    
    return($element_headers);

}

my %options;
my $results = GetOptions( \%options,'inputMatrix|i=s','verbose|v','output|o=s','minDistance|minDist=i','maxDistance|maxDist=i','excludeCis|ec','excludeTrans|et','elementBedFiles|ebf=s@','y_elementBedFiles|yebf=s@','x_elementBedFiles|xebf=s@','outputSuffix|os=s','zoomCoordinates|z=s@','y_zoomCoordinates|yz=s@','x_zoomCoordinates|xz=s@','elementExtension|ee=s') or croak help();
my ($ret,$inputMatrix,$verbose,$output,$minDistance,$maxDistance,$excludeCis,$excludeTrans,$elementBedFiles,$y_elementBedFiles,$x_elementBedFiles,$outputSuffix,$zoomCoordinates,$y_zoomCoordinates,$x_zoomCoordinates,$elementExtension)=check_options( \%options );

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
my ($y_element_headers)=processElements('y',$y_elementBedFiles,$outputSuffix,$elementExtension,$headerBedFile,$verbose) if(@{$y_elementBedFiles} > 0);
my ($x_element_headers)=processElements('x',$x_elementBedFiles,$outputSuffix,$elementExtension,$headerBedFile,$verbose) if(@{$x_elementBedFiles} > 0);

my $element_header2inc={};
my $element_inc2header={};

my $num_element_yHeaders=0;
for(my $y=0;$y<$numYHeaders;$y++) {
    my $yHeader=$inc2header->{ y }->{$y};
    if(exists($y_element_headers->{$yHeader})) {
        $element_header2inc->{ y }->{$yHeader}=$num_element_yHeaders;
        $element_inc2header->{ y }->{$num_element_yHeaders}=$yHeader;
        $num_element_yHeaders++;
    }
}
my $num_element_xHeaders=0;
for(my $x=0;$x<$numXHeaders;$x++) {
    my $xHeader=$inc2header->{ x }->{$x};
    if(exists($x_element_headers->{$xHeader})) {
        $element_header2inc->{ x }->{$xHeader}=$num_element_xHeaders;
        $element_inc2header->{ x }->{$num_element_xHeaders}=$xHeader;
        $num_element_xHeaders++;
    }
    
}
$subset_header2inc->{ y }=$element_header2inc->{ y } if($num_element_yHeaders > 0);
$subset_inc2header->{ y }=$element_inc2header->{ y } if($num_element_yHeaders > 0);
$subset_header2inc->{ x }=$element_header2inc->{ x } if($num_element_xHeaders > 0);
$subset_inc2header->{ x }=$element_inc2header->{ x } if($num_element_xHeaders > 0);

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

my $zoomCoordinate_str=join('__', @$zoomCoordinates);
$zoomCoordinate_str =~ s/\:/-/g;
my $x_zoomCoordinate_str=join('__', @$x_zoomCoordinates);
$x_zoomCoordinate_str =~ s/\:/-/g;
my $y_zoomCoordinate_str=join('__', @$y_zoomCoordinates);
$y_zoomCoordinate_str =~ s/\:/-/g;
if($x_zoomCoordinate_str eq $y_zoomCoordinate_str) {
    $output .= "---".$zoomCoordinate_str if($zoomCoordinate_str ne "");
} else {
    $output .= "--x-".$x_zoomCoordinate_str if($x_zoomCoordinate_str ne "");
    $output .= "--y-".$y_zoomCoordinate_str if($y_zoomCoordinate_str ne "");
}

$output .= "__".$outputSuffix if(($inputMatrixName eq $output) or ($outputSuffix ne ""));

#read Matrix
my $matrix={};
($matrix)=getData($inputMatrix,$matrixObject,$verbose,$minDistance,$maxDistance,$excludeCis,$excludeTrans);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

my $commentLine=getScriptOpts($ret,$tool);

my $subsetMatrixFile=$output.".subset.matrix.gz";
print STDERR "Writing matrix to file ($subsetMatrixFile)...\n" if($verbose);
writeMatrix($matrix,$inc2header,$subsetMatrixFile,"NA",$commentLine);
print STDERR "\tcomplete\n\n" if($verbose);

