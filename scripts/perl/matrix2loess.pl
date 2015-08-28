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

    my ($inputMatrix,$verbose,$output,$loessObjectFile,$excludeCis,$excludeTrans,$cisAlpha,$cisApproximateFactor,$disableIQRFilter,$minDistance,$maxDistance,$excludeZero,$supressMatrixFiles,$supressLoessFile,$debugMode);
    
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
    
    if( exists($opts->{ loessObjectFile }) ) {
        $loessObjectFile = $opts->{ loessObjectFile };
    } else {
        $loessObjectFile="";
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
    
    if( exists($opts->{ supressMatrixFiles }) ) {
        $supressMatrixFiles = 1;
    } else {
        $supressMatrixFiles = 0;
    }
    
    if( exists($opts->{ supressLoessFile }) ) {
        $supressLoessFile = 1;
    } else {
        $supressLoessFile = 0;
    }
    
    if( exists($opts->{ debugMode }) ) {
        $debugMode = 1;
    } else {
        $debugMode = 0;
    }
    
    if(($excludeCis + $excludeTrans) == 2) {
        print STDERR "\nERROR: must choose to exclude CIS data OR TRANS data (NOT BOTH)...\n";
        help();
    }
    
    return($inputMatrix,$verbose,$output,$loessObjectFile,$excludeCis,$excludeTrans,$cisAlpha,$cisApproximateFactor,$disableIQRFilter,$minDistance,$maxDistance,$excludeZero,$supressMatrixFiles,$supressLoessFile,$debugMode);
}

sub intro() {
    print STDERR "\n";
    
    print STDERR "Tool:\t\tmatrix2loess.pl\n";
    print STDERR "Version:\t".$cworld::dekker::VERSION."\n";
    print STDERR "Summary:\tcalculate the loess (expected/stdev/zScore) for a given matrix\n";
    
    print STDERR "\n";
}

sub help() {
    intro();
    
    print STDERR "Usage: perl matrix2loess.pl [OPTIONS] -i <inputMatrix>\n";
    
    print STDERR "\n";
    
    print STDERR "Required:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-i", "[]", "input matrix file");
    
    print STDERR "\n";
    
    print STDERR "Options:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-v", "[]", "FLAG, verbose mode");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-o", "[]", "prefix for output file(s)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--lof", "[]", "optional loess object file (pre-calculated loess)");   
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ec", "[]", "FLAG, exclude CIS data");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--et", "[]", "FLAG, exclude TRANS data");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ez", "[]", "FLAG, ignore zeros from all calculations");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ca", "[0.01]", "lowess alpha value, fraction of datapoints to smooth over");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--caf", "[1000]", "cis approximate factor to speed up loess, genomic distance / -caffs");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--dif", "[]", "FLAG, disable loess IQR (outlier) filter");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--minDist", "[]", "minimum allowed interaction distance, exclude < N distance (in BP)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--maxDist", "[]", "maximum allowed interaction distance, exclude > N distance (in BP)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ez", "[]", "FLAG, ignore 0s in all calculations");   
    printf STDERR ("\t%-10s %-10s %-10s\n", "--smf", "[]", "FLAG, supress all output matrix files");   
    printf STDERR ("\t%-10s %-10s %-10s\n", "--slf", "[]", "FLAG, supress loess file");   
    printf STDERR ("\t%-10s %-10s %-10s\n", "--d", "[]", "FLAG, enable debug mode");   
    
    print STDERR "\n";
    
    print STDERR "Notes:";
    print STDERR "
    This script calculates a lowess smoothing of the data per distance.
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

sub calculateExpected($$$$$;$) {
    my $matrixObject=shift;
    my $loess=shift;
    my $numYHeaders=shift;
    my $numXHeaders=shift;
    my $matrix=shift;
    #optional
    my $cisApproximateFactor=shift;
    
    my $inc2header=$matrixObject->{ inc2header };
    my $missingValue=$matrixObject->{ missingValue };  
    my $verbose=$matrixObject->{ verbose };
    my $output=$matrixObject->{ output };
    
    my $stdevMatrix={};
    my $expectedMatrix={};

    my ($y,$x);
    for(my $y=0;$y<$numYHeaders;$y++) {
        my $yHeader=$inc2header->{ y }->{$y};
        my $yHeaderObject=getHeaderObject($yHeader);
        for(my $x=0;$x<$numXHeaders;$x++) {
            my $xHeader=$inc2header->{ x }->{$x};    
            my $xHeaderObject=getHeaderObject($xHeader);
            
            my $interactionDistance=getInteractionDistance($matrixObject,$yHeaderObject,$xHeaderObject,$cisApproximateFactor);
            my $loessStdev="NA";
            my $loessValue="NA";
            
            my $score=$matrixObject->{ missingValue };
            $score=$matrix->{$y}->{$x} if(defined($matrix->{$y}->{$x}));
            
            if( ($score ne "NA") and (exists($loess->{$interactionDistance}->{"loess"})) ) {
                $loessValue=$loess->{$interactionDistance}->{"loess"} if($loess->{$interactionDistance}->{"loess"} ne "NA");
                $loessStdev=$loess->{$interactionDistance}->{"stdev"} if($loess->{$interactionDistance}->{"stdev"} ne "NA");
            }
            
            $expectedMatrix->{$y}->{$x}=$loessValue if($loessValue ne "NA");
            $stdevMatrix->{$y}->{$x}=$loessStdev if($loessStdev ne "NA");
                                        
        }
        my $pcComplete = round((($y/($numYHeaders-1))*100),2);
        print STDERR "\e[A" if(($verbose) and ($y != 0));
        printf STDERR "\t%.2f%% complete (".$y."/".($numYHeaders-1).")...\n", $pcComplete if($verbose);
    }
    
    my $stdevFile=$output.".stdev.matrix.gz";
    writeMatrix($stdevMatrix,$inc2header,$stdevFile,"NA");

    my $expectedFile=$output.".expected.matrix.gz";
    writeMatrix($expectedMatrix,$inc2header,$expectedFile,"NA");

}

my %options;
my $results = GetOptions( \%options,'inputMatrix|i=s','verbose|v','output|o=s','loessObjectFile|lof=s','excludeCis|ec','excludeTrans|et','cisAlpha|ca=f','cisApproximateFactor|caf=i','disableIQRFilter|dif','minDistance|minDist=i','maxDistance|maxDist=i','excludeZero|ez','supressMatrixFiles|smf','supressLoessFile|slf','debugMode|d') or croak help();

my ($inputMatrix,$verbose,$output,$loessObjectFile,$excludeCis,$excludeTrans,$cisAlpha,$cisApproximateFactor,$disableIQRFilter,$minDistance,$maxDistance,$excludeZero,$supressMatrixFiles,$supressLoessFile,$debugMode)=check_options( \%options );

my $includeCis=flipBool($excludeCis);
my $includeTrans=flipBool($excludeTrans);

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
my $loessFile=$output.$loessMeta.".loess.gz";
$loessObjectFile=$output.$loessMeta.".loess.object.gz" if($loessObjectFile eq "");

my $inputDataCis=[];
my $inputDataTrans=[];
if(!validateLoessObject($loessObjectFile)) {
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
$loess=calculateLoess($matrixObject,$inputDataCis,$loessFile,$loessObjectFile,$cisAlpha,$disableIQRFilter,$excludeZero);

# plot cis-expected
system("Rscript '".$scriptPath."/R/plotLoess.R' '".$cwd."' '".$loessFile."' > /dev/null") if(scalar @{ $inputDataCis } > 0);

# calculate cis-expected
$loess=calculateTransExpected($inputDataTrans,$excludeZero,$loess,$loessObjectFile,$verbose);

#get zScore matrix
print STDERR "calculating z-score matrix...\n" if($verbose);
my $zScoreMatrix=calculateZscore($matrixObject,$matrix,$loess,$cisApproximateFactor,$excludeZero);
my $zScoreFile=$output.".zScore.matrix.gz";
writeMatrix($zScoreMatrix,$inc2header,$zScoreFile,"NA");
undef $zScoreMatrix;
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

# exit here is matrices/plots are not needed
exit if($supressMatrixFiles);

# get expected/stdev matrix
print STDERR "calculating expected matrix...\n" if($verbose);
calculateExpected($matrixObject,$loess,$numYHeaders,$numXHeaders,$matrix,$cisApproximateFactor);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

#get log2ratio matrix
print STDERR "calculating log2ratio matrix...\n" if($verbose);
my $log2ratioMatrix=calculateLog2Ratio($matrixObject,$matrix,$loess,$cisApproximateFactor,$excludeZero);
my $log2ratioFile=$output.".log2ratio.matrix.gz";
writeMatrix($log2ratioMatrix,$inc2header,$log2ratioFile,"NA");
print STDERR "\tdone\n" if($verbose);

#undef log2ratio matrix
undef $log2ratioMatrix;

print STDERR "\n" if($verbose);

#get obs-exp matrix
print STDERR "calculating obsMinusExp matrix...\n" if($verbose);
my $obsMinusExpMatrix=calculateObsMinusExp($matrixObject,$matrix,$loess,$cisApproximateFactor,$excludeZero);
my $obsMinusExpFile=$output.".obs-exp.matrix.gz";
writeMatrix($obsMinusExpMatrix,$inc2header,$obsMinusExpFile,"NA");
print STDERR "\tdone\n" if($verbose);

#undef obsMinusExp matrix
undef $obsMinusExpMatrix;

print STDERR "\n" if($verbose);