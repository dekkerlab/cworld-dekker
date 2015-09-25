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
    
    my ($inputMatrix,$verbose,$output,$loessObjectFile,$includeCis,$minDistance,$maxDistance,$cisAlpha,$disableIQRFilter,$cisApproximateFactor,$includeTrans,$logTransform,$excludeZero,$factorMode,$inputFactorFile,$debugMode);
 
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
    
    if( exists($opts->{ loessObjectFile }) ) {
        $loessObjectFile = $opts->{ loessObjectFile };
    } else {
        $loessObjectFile="";
    }
    
    if( exists($opts->{ includeCis }) ) {
        $includeCis=1;
    } else {
        $includeCis=0;
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
    
    if( exists($opts->{ includeTrans }) ) {
        $includeTrans=1;
    } else {
        $includeTrans=0;
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
        croak "invalid factor mode (zScore,obsExp,zScore+obsExp)" if(($factorMode ne "zScore") and ($factorMode ne "obsExp") and ($factorMode ne "zScore+obsExp"));
    } else {
        $factorMode="zScore";
    }
    
    if( exists($opts->{ inputFactorFile })) {
        $inputFactorFile = $opts->{ inputFactorFile };
    } else {
        print STDERR "\nERROR: Option inputFactorFile|ff is required.\n";
        help();
    }

    if(($includeCis + $includeTrans) != 1) {
        print STDERR "\nERROR: Must choose to use either include CIS (-ic) data XOR TRANS (-it) data.\n";
        help();
    }

    if( exists($opts->{ debugMode }) ) {
        $debugMode = 1;
    } else {
        $debugMode =0;
    }
    
    $ret->{ inputMatrix }=$inputMatrix;
    $ret->{ verbose }=$verbose;
    $ret->{ output }=$output;
    $ret->{ loessObjectFile }=$loessObjectFile;
    $ret->{ includeCis }=$includeCis;
    $ret->{ minDistance }=$minDistance;
    $ret->{ maxDistance }=$maxDistance;
    $ret->{ cisAlpha }=$cisAlpha;
    $ret->{ disableIQRFilter }=$disableIQRFilter;
    $ret->{ cisApproximateFactor }=$cisApproximateFactor;
    $ret->{ includeTrans }=$includeTrans;
    $ret->{ logTransform }=$logTransform;
    $ret->{ excludeZero }=$excludeZero;
    $ret->{ factorMode }=$factorMode;
    $ret->{ inputFactorFile }=$inputFactorFile;
    $ret->{ debugMode }=$debugMode;
    
    return($ret,$inputMatrix,$verbose,$output,$loessObjectFile,$includeCis,$minDistance,$maxDistance,$cisAlpha,$disableIQRFilter,$cisApproximateFactor,$includeTrans,$logTransform,$excludeZero,$factorMode,$inputFactorFile,$debugMode);
}


sub intro() {
    print STDERR "\n";
    
    print STDERR "Tool:\t\t".$tool."\n";
    print STDERR "Version:\t".$cworld::dekker::VERSION."\n";
    print STDERR "Summary:\tapply correction to factor - external factors\n";
    
    print STDERR "\n";
}

sub help() {
    intro();
    
    print STDERR "Usage: perl applyCorrection.pl [OPTIONS] -i <inputMatrix>\n";
    
    print STDERR "\n";
    
    print STDERR "Required:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-i", "[]", "input matrix file");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ff", "[]", "row/col factor file");
   
    print STDERR "\n";
    
    print STDERR "Options:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-v", "[]", "FLAG, verbose mode");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-o", "[]", "prefix for output file(s)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--lof", "[]", "optional loess object file (pre-calculated loess)");   
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ic", "[]", "model CIS data to detect outlier row/cols");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--it", "[]", "model TRANS data to detect outlier row/cols");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--minDist", "[]", "minimum allowed interaction distance, exclude < N distance (in BP)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--maxDist", "[]", "maximum allowed interaction distance, exclude > N distance (in BP)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ca", "[0.01]", "lowess alpha value, fraction of datapoints to smooth over");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--dif", "[]", "FLAG, disable loess IQR (outlier) filter");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--caf", "[1000]", "cis approximate factor to speed up loess, genomic distance / -caffs");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--lt", "[]", "log transform input data into specified base");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ez", "[]", "FLAG, ignore 0s in all calculations");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--fm", "[]", "outlier detection mode - zScore,obsExp");
   
    print STDERR "\n";
    
    print STDERR "Notes:";
    print STDERR "
    This script can apply a set of primer factors to a matrix [normalize].
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
my $results = GetOptions( \%options,'inputMatrix|i=s','verbose|v','output|o=s','loessObjectFile|lof=s','includeCis|ic','minDistance|minDist=i','maxDistance|maxDist=i','cisAlpha|ca=f','disableIQRFilter|dif','cisApproximateFactor|caf=i','includeTrans|it','logTransform|lt=f','excludeZero|ez','factorMode|fm=s','inputFactorFile|ff=s','debugMode|d') or croak help();
my ($ret,$inputMatrix,$verbose,$output,$loessObjectFile,$includeCis,$minDistance,$maxDistance,$cisAlpha,$disableIQRFilter,$cisApproximateFactor,$includeTrans,$logTransform,$excludeZero,$factorMode,$inputFactorFile,$debugMode)=check_options( \%options );

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
$output=$matrixObject->{ output };

my $excludeCis=flipBool($includeCis);
my $excludeTrans=flipBool($includeTrans);

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

my $rowcolData={};

print STDERR "loading primer factors ...\n" if($verbose);
# load the primer factors
if(($inputFactorFile ne "") and (-e($inputFactorFile))) {
    open (my $IFF , "$inputFactorFile") or croak "Could not open file [$inputFactorFile] - $!";
    while(my $line = <$IFF>) {
        chomp($line);
        next if(($line eq "") or ($line =~ m/^#/));
        
        next if(($line =~ m/^\"/) or ($line =~ /^\s*$/));
        
        my @row = split(/\t/,$line);
        #overwrite the factor data
        $rowcolData->{$row[0]}->{ factor } = $row[1];
    }
}
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

print STDERR "correcting matrix ...\n" if($verbose);
# now normalize the matrix
my $normalizedMatrix=correctMatrix($matrixObject,$matrix,$rowcolData,$includeCis,$includeTrans,$factorMode,$logTransform,$loess,$excludeZero,$cisApproximateFactor);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

print STDERR "writing corrected matrix ...\n" if($verbose);
my $normalizedMatrixFile = $output.".corrected.matrix.gz";
writeMatrix($normalizedMatrix,$inc2header,$normalizedMatrixFile,$missingValue,$commentLine);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);
