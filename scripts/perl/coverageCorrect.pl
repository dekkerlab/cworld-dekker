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
    
    my ($inputMatrix,$verbose,$output,$excludeCis,$excludeTrans,$correctionMode,$minDistance,$maxDistance,$cisAlpha,$disableIQRFilter,$cisApproximateFactor,$logTransform,$excludeZero,$factorMode,$maxIterations,$convergenceThreshold,$debug);
 
    my $ret={};
    
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

    if( exists($opts->{ correctionMode }) ) {
        $correctionMode = $opts->{ correctionMode };
        croak "invalid correction mode [$correctionMode] (cis,trans)" if(($correctionMode ne "cis") and ($correctionMode ne "trans"));
    } else {
        $correctionMode = "trans"
    }

    
    if( exists($opts->{ maxDistance }) ) {
        $maxDistance = $opts->{ maxDistance };
    } else {
        $maxDistance=undef;
    }
    
    if( exists($opts->{ minDistance }) ) {
        $minDistance = $opts->{ minDistance };
    } else {
        $minDistance=undef;
    }
    
    if( exists($opts->{ cisAlpha }) ) {
        $cisAlpha = $opts->{ cisAlpha };
    } else {
        $cisAlpha=0.01;
    }    
    
    if( exists($opts->{ disableIQRFilter }) ) {
        $disableIQRFilter = 1;
    } else {
        $disableIQRFilter = 0;
    }
    
    if( exists($opts->{ cisApproximateFactor }) ) {
        $cisApproximateFactor = $opts->{ cisApproximateFactor };
    } else {
        $cisApproximateFactor=1000;
    }
    
    if( exists($opts->{ logTransform }) ) {
        $logTransform = $opts->{ logTransform };
    } else {
        $logTransform=0;
    }    
    
    if( exists($opts->{ excludeZero }) ) {
        $excludeZero = 1;
    } else {
        $excludeZero = 0;
    }
    
    if( exists($opts->{ factorMode }) ) {
        $factorMode = $opts->{ factorMode };
        croak "invalid factor mode [$factorMode] (zScore,obsExp,zScore+obsExp)" if(($factorMode ne "zScore") and ($factorMode ne "obsExp") and ($factorMode ne "zScore+obsExp"));
    } else {
        $factorMode="zScore";
    }
        
    if( exists($opts->{ maxIterations }) ) {
        $maxIterations = $opts->{ maxIterations };
    } else {
        $maxIterations = 100;
    }
    
    if( exists($opts->{ convergenceThreshold }) ) {
        $convergenceThreshold = $opts->{ convergenceThreshold };
    } else {
        $convergenceThreshold = 0.1;
    }

    if( exists($opts->{ debug }) ) {
        $debug = 1;
    } else {
        $debug =0;
    }
    
    $ret->{ inputMatrix }=$inputMatrix;
    $ret->{ verbose }=$verbose;
    $ret->{ output }=$output;
    $ret->{ excludeCis }=$excludeCis;
    $ret->{ excludeTrans }=$excludeTrans;
    $ret->{ correctionMode }=$correctionMode;
    $ret->{ minDistance }=$minDistance;
    $ret->{ maxDistance }=$maxDistance;
    $ret->{ cisAlpha }=$cisAlpha;
    $ret->{ disableIQRFilter }=$disableIQRFilter;
    $ret->{ cisApproximateFactor }=$cisApproximateFactor;
    $ret->{ logTransform }=$logTransform;
    $ret->{ excludeZero }=$excludeZero;
    $ret->{ factorMode }=$factorMode;
    $ret->{ maxIterations }=$maxIterations;
    $ret->{ convergenceThreshold }=$convergenceThreshold;
    $ret->{ debug }=$debug;
    
    return($ret,$inputMatrix,$verbose,$output,$excludeCis,$excludeTrans,$correctionMode,$minDistance,$maxDistance,$cisAlpha,$disableIQRFilter,$cisApproximateFactor,$logTransform,$excludeZero,$factorMode,$maxIterations,$convergenceThreshold,$debug);
}

sub intro() {
    print STDERR "\n";
    
    print STDERR "Tool:\t\t".$tool."\n";
    print STDERR "Version:\t".$cworld::dekker::VERSION."\n";
    print STDERR "Summary:\tcan perform coverage correction on matrix [balancing] cis/trans\n";
    
    print STDERR "\n";
}

sub help() {
    intro();
    
    print STDERR "Usage: perl coverageCorrect.pl [OPTIONS] -i <inputMatrix>\n";
    
    print STDERR "\n";
    
    print STDERR "Required:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-i", "[]", "input matrix file");
    
    print STDERR "\n";
    
    print STDERR "Options:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-v", "[]", "FLAG, verbose mode");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-o", "[]", "prefix for output file(s)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--cm", "[]", "correction mode (cis,trans), use cis or trans data for bin/primer performance");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ec", "[]", "FLAG, exclude CIS data");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--et", "[]", "FLAG, exclude TRANS data");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--minDist", "[]", "minimum allowed interaction distance, exclude < N distance (in BP)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--maxDist", "[]", "maximum allowed interaction distance, exclude > N distance (in BP)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ca", "[0.01]", "lowess alpha value, fraction of datapoints to smooth over");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--dif", "[]", "FLAG, disable loess IQR (outlier) filter");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--caf", "[1000]", "cis approximate factor to speed up loess, genomic distance / -caffs");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--lt", "[]", "log transform input data into specified base");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ez", "[]", "FLAG, ignore 0s in all calculations");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--fm", "[zScore]", "factorMode, how to apply the correction factors (zScore,obsExp)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--mi", "[]", "max iterations allowed");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ct", "[0.1]", "optional, convergance threshold");
    print STDERR "\n";
    
    print STDERR "Notes:";
    print STDERR "
    This script can compare two matrices via scatter plot regression.
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

sub getDeltaFactors($$) {
    my $rowcolData=shift;
    my $lastPrimerData=shift;

    my @tmpDeltaArr=();
    my @tmpFactorArr=();
    foreach my $primer (keys %$rowcolData) {
        my $factor=$rowcolData->{$primer}->{ factor };
        my $lastFactor=$lastPrimerData->{$primer}->{ factor };
    
        next if(($factor eq "NA") or ($lastFactor eq "NA"));

        my $delta=abs($factor-$lastFactor);
        push(@tmpDeltaArr,$delta);
        push(@tmpFactorArr,abs($factor));
    }
    my $tmpDeltaArrStats=listStats(\@tmpDeltaArr);
    my $maxDelta=$tmpDeltaArrStats->{ max };
    my $meanDelta=$tmpDeltaArrStats->{ iqrMean };
    
    my $tmpFactorArrStats=listStats(\@tmpFactorArr);
    my $maxFactor=$tmpFactorArrStats->{ max };
    my $meanFactor=$tmpFactorArrStats->{ iqrMean };
    
    return($meanDelta);
    
}


my %options;
my $results = GetOptions( \%options,'inputMatrix|i=s','verbose|v','output|o=s','excludeCis|ec','excludeTrans|et','correctionMode|cm=s','minDistance|minDist=i','maxDistance|maxDist=i','cisAlpha|ca=f','disableIQRFilter|dif','cisApproximateFactor|caf=i','logTransform|lt=f','excludeZero|ez','factorMode|fm=s','maxIterations|mi=i','convergenceThreshold|ct=f','debug|d') or croak help();
my ($ret,$inputMatrix,$verbose,$output,$excludeCis,$excludeTrans,$correctionMode,$minDistance,$maxDistance,$cisAlpha,$disableIQRFilter,$cisApproximateFactor,$logTransform,$excludeZero,$factorMode,$maxIterations,$convergenceThreshold,$debug)=check_options( \%options );

intro() if($verbose);

#get the absolute path of the script being used
my $cwd = getcwd();
my $fullScriptPath=abs_path($0);
my @fullScriptPathArr=split(/\//,$fullScriptPath);
@fullScriptPathArr=@fullScriptPathArr[0..@fullScriptPathArr-3];
my $scriptPath=join("/",@fullScriptPathArr);
my $commentLine=getScriptOpts($ret,$tool);

croak "inputMatrix [$inputMatrix] does not exist" if(!(-e $inputMatrix));

# get matrix information
my $matrixObject=getMatrixObject($inputMatrix,$output,$verbose);
my $inc2header=$matrixObject->{ inc2header };
my $header2inc=$matrixObject->{ header2inc };
my $numYHeaders=$matrixObject->{ numYHeaders };
my $numXHeaders=$matrixObject->{ numXHeaders };
my $missingValue=$matrixObject->{ missingValue };
my $symmetrical=$matrixObject->{ symmetrical };
my $inputMatrixName=$matrixObject->{ inputMatrixName };

my $includeCis=flipBool($excludeCis);
my $includeTrans=flipBool($excludeTrans);

my $factor_includeCis=0;
$factor_includeCis=1 if($correctionMode eq "cis");
my $factor_includeTrans=0;
$factor_includeTrans=1 if($correctionMode eq "trans");

my $tmpInputMatrix = "input__".$inputMatrixName.".matrix.gz";
$tmpInputMatrix=$output if($output ne "");
$output=$matrixObject->{ output };

my $tmpInputMatrixName = getFileName($tmpInputMatrix);
system("cp '".$inputMatrix."' '".$tmpInputMatrix."'");

my $rowcolData={};
my $lastPrimerData={};
my $iterationKey="";

my $iterationNumber=0;
my $meanDelta=1;

my $iterationFile=$output.".log";
open(LOG,outputWrapper($iterationFile)) or croak "Could not open file [$iterationFile] - $!";

print STDERR "\n";

while(($meanDelta > $convergenceThreshold) and ($iterationNumber <= $maxIterations)) {
    
    $iterationKey="i".$iterationNumber."__" if($debug);
    $tmpInputMatrixName = $iterationKey.getFileName($tmpInputMatrix) if($iterationNumber > 0);
    
    print STDERR "correction iteration # $iterationNumber ($tmpInputMatrixName) ...\n";
    
    # get matrix information
    my $matrixObject=getMatrixObject($tmpInputMatrix,$output,$verbose);
    my $inc2header=$matrixObject->{ inc2header };
    my $header2inc=$matrixObject->{ header2inc };
    my $numYHeaders=$matrixObject->{ numYHeaders };
    my $numXHeaders=$matrixObject->{ numXHeaders };
    my $missingValue=$matrixObject->{ missingValue };

    # get matrix data
    my ($matrix)=getData($tmpInputMatrix,$matrixObject,$verbose,$minDistance,$maxDistance,$excludeCis,$excludeTrans);

    print STDERR "\n" if($verbose);

    # calculate LOWESS smoothing for cis data

     my $loessMeta="";
    $loessMeta .= "--ic" if($includeCis);
    $loessMeta .= "--it" if($includeTrans);
    $loessMeta .= "--maxDist".$maxDistance if(defined($maxDistance));
    $loessMeta .= "--ez" if($excludeZero);
    $loessMeta .= "--caf".$cisApproximateFactor;
    $loessMeta .= "--ca".$cisAlpha;
    $loessMeta .= "--dif" if($disableIQRFilter);
    my $loessFile=$tmpInputMatrixName.$loessMeta.".loess.gz";
    my $loessObjectFile=$tmpInputMatrixName.$loessMeta.".loess.object.gz";
    system("rm '".$loessObjectFile."'") if(-e($loessObjectFile));
    
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
    
    print STDERR "calculating y-axis factors ... \n" if($verbose);
    my $yFactorFile=$tmpInputMatrixName.".y.factors";
    my $normalPrimerData=getRowColFactor($matrixObject,$tmpInputMatrix,$tmpInputMatrixName,$factor_includeCis,$minDistance,$maxDistance,$factor_includeTrans,$logTransform,$loess,$factorMode,$yFactorFile,$inputMatrixName,$cisApproximateFactor,$excludeZero);
   
    print STDERR "\n" if($verbose);

    print STDERR "transposing matrix ... \n" if($verbose);
    my $transposedMatrix=transposeMatrix($tmpInputMatrix);
    my $transposedMatrixObject=getMatrixObject($transposedMatrix);
    
    print STDERR "\n" if($verbose);
    
    print STDERR "calculating x-axis factors ... \n" if($verbose);
    my $xFactorFile=$tmpInputMatrixName.".x.factors";
    my $transposedPrimerData=getRowColFactor($transposedMatrixObject,$transposedMatrix,$tmpInputMatrixName,$factor_includeCis,$minDistance,$maxDistance,$factor_includeTrans,$logTransform,$loess,$factorMode,$xFactorFile,$inputMatrixName,$cisApproximateFactor,$excludeZero);
    
    print STDERR "\n" if($verbose);

    # cleanup
    system("rm ".$transposedMatrix);
    system("rm ".$yFactorFile);
    system("rm ".$xFactorFile);
    
    # validate primer zScore hashes
    foreach my $yPrimer (keys %$normalPrimerData) {
        my $yPrimerFactor=$normalPrimerData->{$yPrimer}->{ factor };
        next if($yPrimerFactor eq "NA");
        if(exists($transposedPrimerData->{$yPrimer})) {
            croak "matrix headers symmetrical, but data is not ($yPrimer)" if($transposedPrimerData->{$yPrimer}->{ factor } != $yPrimerFactor);
        }
    }
    foreach my $xPrimer (keys %$transposedPrimerData) {
        my $xPrimerFactor=$transposedPrimerData->{$xPrimer}->{ factor };
        next if($xPrimerFactor eq "NA");
        if(exists($normalPrimerData->{$xPrimer})) {
            croak "matrix headers symmetrical, but data is not ($xPrimer)" if($normalPrimerData->{$xPrimer}->{ factor } != $xPrimerFactor);
        }
    }

    # now combine x and y hashes - symmetrical is AOK
    my %primerDataHash=();
    %primerDataHash=(%$normalPrimerData,%$transposedPrimerData);
    $rowcolData=\%primerDataHash;

    # calculate and draw histogram 
    my $allFactorFile=$tmpInputMatrixName.".mean-zScores";

    # re-scale and output
    open(OUT,outputWrapper($allFactorFile)) or croak "Could not open file [$allFactorFile] - $!";
    foreach my $primerName ( keys %$rowcolData ) {
        my $factor=$rowcolData->{$primerName}->{ factor };
        next if($factor eq "NA");
        print OUT "$primerName\t$factor\n";
    }
    close(OUT);

    # plot mean-zsocres
    system("Rscript '".$scriptPath."/R/factorHistogram.R' '".$cwd."' '".$allFactorFile."' '".$tmpInputMatrixName."' > /dev/null");

    # now normalize the matrix
    my $normalizedMatrix=correctMatrix($matrixObject,$matrix,$rowcolData,$includeCis,$includeTrans,$factorMode,$logTransform,$loess,$excludeZero,$cisApproximateFactor);
    
    my $normalizedMatrixFile = $iterationKey.$output.".corrected.matrix.gz";
    writeMatrix($normalizedMatrix,$inc2header,$normalizedMatrixFile,$missingValue,$commentLine);
    
    $meanDelta=getDeltaFactors($rowcolData,$lastPrimerData) if($iterationNumber != 0);
    print STDERR "\tdelta = $meanDelta\n";
    print STDERR "\n" if($verbose);
    print LOG "iteration #".$iterationNumber."\t$meanDelta\t$convergenceThreshold\n";
    
    $tmpInputMatrix=$normalizedMatrixFile;
    $tmpInputMatrixName=getFileName($tmpInputMatrix);
    
    $lastPrimerData=$rowcolData;
    $iterationNumber++;
}

print STDERR "\nconverged in $iterationNumber iterations!\n";
