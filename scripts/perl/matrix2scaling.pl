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

    my ($inputMatrixArray,$verbose,$output,$loessObjectFile,$cisAlpha,$disableIQRFilter,$minDistance,$maxDistance,$cisApproximateFactor,$excludeZero);
    
    my $ret={};
    
    if( exists($opts->{ inputMatrixArray }) ) {
        $inputMatrixArray = $opts->{ inputMatrixArray };
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
    
    if( exists($opts->{ loessObjectFile }) ) {
        $loessObjectFile = $opts->{ loessObjectFile };
    } else {
        $loessObjectFile="";
    }
    
    if( exists($opts->{ cisAlpha }) ) {
        $cisAlpha = $opts->{ cisAlpha };
    } else {        
        $cisAlpha=0.01;
    }    
    
    if( exists($opts->{ cisApproximateFactor }) ) {
        $cisApproximateFactor = $opts->{ cisApproximateFactor };
    } else {
        $cisApproximateFactor=1000;
    }
    
    
    if( exists($opts->{ disableIQRFilter }) ) {
        $disableIQRFilter = 1;
    } else {
        $disableIQRFilter=0;
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
    
    if( exists($opts->{ excludeZero }) ) {
        $excludeZero = 1;
    } else {
        $excludeZero = 0;
    }
    
    $ret->{ inputMatrixArray }=$inputMatrixArray;
    $ret->{ verbose }=$verbose;
    $ret->{ output }=$output;
    $ret->{ loessObjectFile }=$loessObjectFile;
    $ret->{ cisAlpha }=$cisAlpha;
    $ret->{ disableIQRFilter }=$disableIQRFilter;
    $ret->{ minDistance }=$minDistance;
    $ret->{ maxDistance }=$maxDistance;
    $ret->{ cisApproximateFactor }=$cisApproximateFactor;
    $ret->{ excludeZero }=$excludeZero;
    
    return($ret,$inputMatrixArray,$verbose,$output,$loessObjectFile,$cisAlpha,$disableIQRFilter,$minDistance,$maxDistance,$cisApproximateFactor,$excludeZero);
}

sub intro() {
    print STDERR "\n";
    
    print STDERR "Tool:\t\t".$tool."\n";
    print STDERR "Version:\t".$cworld::dekker::VERSION."\n";
    print STDERR "Summary:\ttransform matrix into scaling (polymer) plot\n";
    
    print STDERR "\n";
}

sub help() {
    intro();
    
    print STDERR "Usage: perl matrix2scaling.pl [OPTIONS] -i <inputMatrix>\n";
    
    print STDERR "\n";
    
    print STDERR "Required:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-i", "[]", "input matrix file (can accept multiple files)");
    
    print STDERR "\n";
    
    print STDERR "Options:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-v", "[]", "FLAG, verbose mode");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-o", "[]", "prefix for output file(s)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--lof", "[]", "optional loess object file (pre-calculated loess)");   
    printf STDERR ("\t%-10s %-10s %-10s\n", "--n", "[]", "optional job name for output files");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ca", "[0.01]", "lowess alpha value, fraction of datapoints to smooth over");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--dif", "[]", "FLAG, disable loess IQR (outlier) filter");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--minDist", "[]", "minimum allowed interaction distance, exclude < N distance (in BP)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--maxDist", "[]", "maximum allowed interaction distance, exclude > N distance (in BP)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--caf", "[1000]", "cis approximate factor to speed up loess, genomic distance / -caffs");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ez", "[]", "FLAG, ignore 0s in all calculations");
    
    print STDERR "\n";
    
    print STDERR "Notes:";
    print STDERR "
    This script calculates a scaling plot of a given matrix.
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

my %options;
my $results = GetOptions( \%options,'inputMatrixArray|i=s@','verbose|v','output|o=s','loessObjectFile|lof=s','cisAlpha|ca=f','disableIQRFilter|dif=s','minDistance|minDist=i','maxDistance|maxDist=i','cisApproximateFactor|caf=i','excludeZero|ez') or croak help();
my ($ret,$inputMatrixArray,$verbose,$output,$loessObjectFile,$cisAlpha,$disableIQRFilter,$minDistance,$maxDistance,$cisApproximateFactor,$excludeZero)=check_options( \%options );

intro() if($verbose);

#get the absolute path of the script being used
my $cwd = getcwd();
my $fullScriptPath=abs_path($0);
my @fullScriptPathArr=split(/\//,$fullScriptPath);
@fullScriptPathArr=@fullScriptPathArr[0..@fullScriptPathArr-3];
my $scriptPath=join("/",@fullScriptPathArr);
my $commentLine=getScriptOpts($ret,$tool);

my $nFiles=@{$inputMatrixArray};
$output="myScaling" if(($nFiles > 1) and ($output eq ""));
$output=getFileName($inputMatrixArray->[0]) if(($nFiles) and ($output eq ""));

# global options
my $excludeCis=0;
my $excludeTrans=0;
my $includeCis=flipBool($excludeCis);
my $includeTrans=flipBool($excludeTrans);

my $filesToPlot="";
my $outputsToPlot="";
my $averageTransSignals="";

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
    my $tmp_output=$inputMatrixName;
        
    # get matrix data
    my ($matrix)=getData($inputMatrix,$matrixObject,$verbose,$minDistance,$maxDistance,$excludeCis,$excludeTrans);

    print STDERR "\n" if($verbose);

    # calculate LOWESS smoothing for cis data

    croak "loessObjectFile [$loessObjectFile] does not exist" if( ($loessObjectFile ne "") and (!(-e $loessObjectFile)) );

    my $loessMeta="";
    $loessMeta .= "--ic" if($includeCis);
    $loessMeta .= "--it" if($includeTrans);
    $loessMeta .= "--maxDist".$maxDistance if(defined($maxDistance));
    $loessMeta .= "--ez" if($excludeZero);
    $loessMeta .= "--caf".$cisApproximateFactor;
    $loessMeta .= "--ca".$cisAlpha;
    $loessMeta .= "--dif" if($disableIQRFilter);
    my $tmp_loessFile=$tmp_output.$loessMeta.".loess.gz";
    my $tmp_loessObjectFile=$tmp_output.$loessMeta.".loess.object.gz";
    $tmp_loessObjectFile=$loessObjectFile if($loessObjectFile ne "");
    
    my $inputDataCis=[];
    my $inputDataTrans=[];
    if(!validateLoessObject($tmp_loessObjectFile)) {
        # dump matrix data into input lists (CIS + TRANS)
        print STDERR "seperating cis/trans data...\n" if($verbose);
        ($inputDataCis,$inputDataTrans)=matrix2inputlist($matrixObject,$matrix,$includeCis,$includeTrans,$minDistance,$maxDistance,$excludeZero,$cisApproximateFactor);
        croak "$inputMatrixName - no avaible CIS data" if((scalar @{ $inputDataCis } <= 0) and ($includeCis) and ($includeTrans == 0));
        croak "$inputMatrixName - no avaible TRANS data" if((scalar @{ $inputDataTrans } <= 0) and ($includeTrans) and ($includeCis == 0));
        print STDERR "\n" if($verbose);
    }

    # init loess object [hash]
    my $loess={};

    # calculate cis-expected
    $loess=calculateLoess($matrixObject,$inputDataCis,$tmp_loessFile,$tmp_loessObjectFile,$cisAlpha,$disableIQRFilter,$excludeZero);

    # plot cis-expected
    system("Rscript '".$scriptPath."/R/plotLoess.R' '".$cwd."' '".$tmp_loessFile."' > /dev/null") if(scalar @{ $inputDataCis } > 0);

    # calculate cis-expected
    $loess=calculateTransExpected($inputDataTrans,$excludeZero,$loess,$tmp_loessObjectFile,$verbose);

    my $averageTrans="NA";
    $averageTrans=$loess->{-1}->{ loess } if(exists($loess->{-1}->{ loess }));
    
    $filesToPlot .= ",".$tmp_loessFile if($filesToPlot ne "");
    $filesToPlot = $tmp_loessFile if($filesToPlot eq "");
    
    $outputsToPlot .= ",".getFileName($tmp_loessFile) if($outputsToPlot ne "");
    $outputsToPlot = getFileName($tmp_loessFile) if($outputsToPlot eq "");
    
    $averageTransSignals .= ",".$averageTrans if($averageTransSignals ne "");
    $averageTransSignals = $averageTrans if($averageTransSignals eq "");

}

my @filesToPlotArray=split(/,/,$filesToPlot);
my @averageTransSignalsArray=split(/,/,$averageTransSignals);

croak "cis/trans mismatch" if(@filesToPlotArray != @averageTransSignalsArray);

my $scalingPlotScriptPath=$scriptPath."/R/plotScaling.R";
system("Rscript '".$scalingPlotScriptPath."' '".$cwd."' '".$filesToPlot."' '".$outputsToPlot."' '".$averageTransSignals."' '".$output."' &> /dev/null");

# do cleanup
#for(my $i=0;$i<@filesToPlotArray;$i++) {
#    my $file=$filesToPlotArray[$i];
#   system("rm '$file'");
#}