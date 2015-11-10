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
    
    my ($inputMatrix,$verbose,$output,$loessObjectFile,$manualOutlierFile,$minDistance,$maxDistance,$cisAlpha,$disableIQRFilter,$cisApproximateFactor,$factorMode,$includeCis,$includeTrans,$iqrMultiplier,$logTransform,$excludeZero);
    
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
    
    if( exists($opts->{ manualOutlierFile }) ) {
        $manualOutlierFile = $opts->{ manualOutlierFile };
    } else {
        $manualOutlierFile="";
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
    
    if( exists($opts->{ factorMode }) ) {
        $factorMode = $opts->{ factorMode };
        croak "invalid factor mode (zScore,obsExp)" if(($factorMode ne "zScore") and ($factorMode ne "obsExp"));
    } else {
        $factorMode="obsExp";
    }

    if( exists($opts->{ includeCis }) ) {
        $includeCis=1;
    } else {
        $includeCis=0;
    }
    
    if( exists($opts->{ includeTrans }) ) {
        $includeTrans=1;
    } else {
        $includeTrans=0;
    }
    
    if( exists($opts->{ iqrMultiplier }) ) {
        $iqrMultiplier = $opts->{ iqrMultiplier };
    } else {
        $iqrMultiplier = 1.5;
    }    
    
    if( exists($opts->{ logTransform }) ) {
        $logTransform = $opts->{ logTransform };
    } else {
        $logTransform = 0;
    }    
    
    if( exists($opts->{ excludeZero }) ) {
        $excludeZero = 1;
    } else {
        $excludeZero = 0;
    }
    
    if(($manualOutlierFile eq "") and (($includeCis + $includeTrans) != 1)) {
        print STDERR "\nERROR: must choose either to include CIS (-ic) data OR TRANS (-it) data.\n";
        help();
    }
    
    $ret->{ inputMatrix }=$inputMatrix;
    $ret->{ verbose }=$verbose;
    $ret->{ output }=$output;
    $ret->{ loessObjectFile }=$loessObjectFile;
    $ret->{ manualOutlierFile }=$manualOutlierFile;
    $ret->{ minDistance }=$minDistance;
    $ret->{ maxDistance }=$maxDistance;
    $ret->{ cisAlpha }=$cisAlpha;
    $ret->{ disableIQRFilter }=$disableIQRFilter;
    $ret->{ cisApproximateFactor }=$cisApproximateFactor;
    $ret->{ factorMode }=$factorMode;
    $ret->{ includeCis }=$includeCis;
    $ret->{ includeTrans }=$includeTrans;
    $ret->{ iqrMultiplier }=$iqrMultiplier;
    $ret->{ logTransform }=$logTransform;
    $ret->{ excludeZero }=$excludeZero;
    
    return($ret,$inputMatrix,$verbose,$output,$loessObjectFile,$manualOutlierFile,$minDistance,$maxDistance,$cisAlpha,$disableIQRFilter,$cisApproximateFactor,$factorMode,$includeCis,$includeTrans,$iqrMultiplier,$logTransform,$excludeZero);
}

sub intro() {
    print STDERR "\n";
    
    print STDERR "Tool:\t\t".$tool."\n";
    print STDERR "Version:\t".$cworld::dekker::VERSION."\n";
    print STDERR "Summary:\tfilters out row/col from C data matrix file\n";
    
    print STDERR "\n";
}

sub help() {
    intro();
    
    print STDERR "Usage: perl anchorPurge.pl [OPTIONS] -i <inputMatrix>\n";
    
    print STDERR "\n";
    
    print STDERR "Required:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-i", "[]", "input matrix file");
    
    print STDERR "\n";
    
    print STDERR "Options:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-v", "[]", "FLAG, verbose mode");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-o", "[]", "prefix for output file(s)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--lof", "[]", "optional loess object file (pre-calculated loess)");   
    printf STDERR ("\t%-10s %-10s %-10s\n", "--im", "[1.5]", "iqr multiplier for setting outlier bounds, q1-<N>*IQR & q3+<N>*IQR");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ic", "[]", "model CIS data to detect outlier row/cols");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--it", "[]", "model TRANS data to detect outlier row/cols");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--mof", "[]", "optional manual outlier file, 1 header per line to be filtered.");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--minDist", "[]", "minimum allowed interaction distance, exclude < N distance (in BP)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--maxDist", "[]", "maximum allowed interaction distance, exclude > N distance (in BP)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ca", "[0.01]", "lowess alpha value, fraction of datapoints to smooth over");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--dif", "[]", "FLAG, disable loess IQR (outlier) filter");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--caf", "[1000]", "cis approximate factor to speed up loess, genomic distance / -caffs");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--lt", "[]", "log transform input data into specified base");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ez", "[]", "FLAG, ignore 0s in all calculations");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--fm", "[obsExp]", "outlier detection mode - (zScore,obsExp)");
    
    print STDERR "\n";
    
    print STDERR "Notes:";
    print STDERR "
    This script detects and removes outlier row/cols.
    Outlier row/cols can be detected by either cis or trans data.
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

sub processManualOutliers($) {
    my $manualOutlierFile=shift;
    
    my %manualOutliers=();
    
    open(IN,inputWrapper($manualOutlierFile)) or croak "Could not open file [$manualOutlierFile] - $!";
    while(my $line = <IN>) {
        chomp($line);
        next if(($line eq "") or ($line =~ m/^#/));
        
        $manualOutliers{$line}=1;
    }
    
    return(\%manualOutliers);
}
    
sub writeOutliers($$$$$$) {
    my $matrixObject=shift;
    my $rowcolData=shift;
    my $output=shift;
    my $botCutoff=shift;
    my $topCutoff=shift;
    my $commentLine=shift;
    
    my $verbose=$matrixObject->{ verbose };
    
    my $toRemoveFile=$output.".toRemove";

    open(OUT,outputWrapper($toRemoveFile,$commentLine)) or croak "Could not open file [$toRemoveFile] - $!";

    my %anchorOutliers=();
    my $nAnchorOutliers=0;
    
    foreach my $primerName ( keys %$rowcolData ) {
        my $meanFactor=$rowcolData->{$primerName}->{ factor };
        next if($meanFactor eq "NA");
        
        if(($meanFactor <= $botCutoff) or ($meanFactor >= $topCutoff)) {
            # log the outliers
            $anchorOutliers{$primerName}=$meanFactor;
            $nAnchorOutliers++;
            print OUT "$primerName\t$meanFactor\n";
            print STDERR "\tremoving $primerName (".$meanFactor.")\n" if($verbose);
        }
    }

    close(OUT);
   
    print STDERR "\tremoved $nAnchorOutliers anchor outliers\n" if($verbose);
    
    return(\%anchorOutliers);
    
}

sub writePrimerFactors($$$$) {
    my $rowcolData=shift;
    my $output=shift;
    my $scriptPath=shift;
    my $commentLine=shift;
    
    # calculate and draw histogram 
    my $allFactorFile=$output.".factor";

    # re-scale and output
    open(FACTOR,outputWrapper($allFactorFile,$commentLine)) or croak "Could not open file [$allFactorFile] - $!";

    foreach my $primerName ( keys %$rowcolData ) {
        my $mean=$rowcolData->{$primerName}->{ factor };
        
        next if($mean eq "NA");
        
        print FACTOR "$primerName\t$mean\n";
    }
    close(FACTOR);
        
    my $cwd = getcwd();
    
    # plot re-scaled values
    system("Rscript '".$scriptPath."/R/factorHistogram.R' '".$cwd."' '".$allFactorFile."' '".$output."' > /dev/null");
    system("rm $allFactorFile");

}

sub calculateBounds($$$$$$) {
    my $matrixObject=shift;
    my $rowcolData=shift;
    my $iqrMultiplier=shift;
    my $output=shift;
    my $scriptPath=shift;
    my $commentLine=shift;
    
    my $verbose=$matrixObject->{ verbose };
    
    # get all the means
    my @tmpAll=();
    my @tmpPos=();
    my @tmpNeg=();
    my $allPrimerFactorFile=$output.".allPrimerFactors";
    open(OUTALL,outputWrapper($allPrimerFactorFile,$commentLine)) or croak "Could not open file [$allPrimerFactorFile] - $!";
    foreach my $primerName ( keys %$rowcolData ) {
        my $mean=$rowcolData->{$primerName}->{ factor };
        next if($mean eq "NA");
        
        push(@tmpAll,$mean);
        print OUTALL "$primerName\t$mean\n";
        
        push(@tmpPos,$mean) if($mean > 0);
        push(@tmpNeg,$mean) if($mean < 0);
    }
    close(OUTALL);

    my $cwd = getcwd();
    
    # plot pos values
    system("Rscript '".$scriptPath."/R/factorHistogram.R' '".$cwd."' '".$allPrimerFactorFile."' '".$output."' > /dev/null");
    
    # calculate bounds
    my $tmpPosArrStats=listStats(\@tmpPos,0.1);
    my $tmpPosBound=$tmpPosArrStats->{ max };
    my $tmpPosMax=$tmpPosArrStats->{ max };
    my $tmpPosMean=$tmpPosArrStats->{ mean };
    my $tmpPosMedian=$tmpPosArrStats->{ median };
    my $tmpPosMin=$tmpPosArrStats->{ min };
    my $tmpPosStdev=$tmpPosArrStats->{ stdev };
    my $tmpPosQ1=$tmpPosArrStats->{ q1 };
    my $tmpPosQ3=$tmpPosArrStats->{ q3 };
    my $tmpPosIQR=$tmpPosArrStats->{ iqr };
    my $tmpPosMAD=$tmpPosArrStats->{ mad };

    my $tmpNegArrStats=listStats(\@tmpNeg,0.1);
    my $tmpNegBound=$tmpNegArrStats->{ min };
    my $tmpNegMax=$tmpNegArrStats->{ max };
    my $tmpNegMean=$tmpNegArrStats->{ mean };
    my $tmpNegMedian=$tmpNegArrStats->{ median };
    my $tmpNegMin=$tmpNegArrStats->{ min };
    my $tmpNegStdev=$tmpNegArrStats->{ stdev };
    my $tmpNegQ1=$tmpNegArrStats->{ q1 };
    my $tmpNegQ3=$tmpNegArrStats->{ q3 };
    my $tmpNegIQR=$tmpNegArrStats->{ iqr };
    my $tmpNegMAD=$tmpNegArrStats->{ mad };
    print STDERR "\tneg\t$tmpNegMin\t$tmpNegQ1\t$tmpNegMedian\t$tmpNegQ3\t$tmpNegMax\t($tmpNegIQR)\t$tmpNegMean\t$tmpNegStdev\n" if($verbose);
    print STDERR "\tpos\t$tmpPosMin\t$tmpPosQ1\t$tmpPosMedian\t$tmpPosQ3\t$tmpPosMax\t($tmpPosIQR)\t$tmpPosMean\t$tmpPosStdev\n" if($verbose);
    
    $tmpPosBound=$tmpPosQ3+($iqrMultiplier*$tmpPosIQR);
    $tmpNegBound=$tmpNegQ1-($iqrMultiplier*$tmpNegIQR);
    print STDERR "\tposBound\t$tmpPosBound\n" if($verbose);
    print STDERR "\tnegBound\t$tmpNegBound\n" if($verbose);
    print STDERR "\tposMax\t$tmpPosMax\n" if($verbose);
    print STDERR "\tnegMin\t$tmpNegMin\n" if($verbose);
    $tmpPosBound = $tmpPosMax if($tmpPosBound > $tmpPosMax); 
    $tmpNegBound = $tmpNegMin if($tmpNegBound < $tmpNegMin); #override to min if outside of bounds
    print STDERR "\tposBound\t$tmpPosBound\n" if($verbose);
    print STDERR "\tnegBound\t$tmpNegBound\n" if($verbose);
    
    my $scaledZScoreLogFile=$output.".re-scaled-zScores.txt";
    open(OUT,outputWrapper($scaledZScoreLogFile,$commentLine)) or croak "Could not open file [$scaledZScoreLogFile] - $!";

    print OUT "zScore bounds\n";
    print OUT "\tposMax\t$tmpPosMax\n";
    print OUT "\tnegMin\t$tmpNegMin\n";
    print OUT "\n";
    print OUT "type\tmin\tmax\tmedian\tmean\tstdev\tq1\tq3\tiqr\n";
    print OUT "\tpos\t$tmpPosMin\t$tmpPosMax\t$tmpPosMedian\t$tmpPosMean\t$tmpPosStdev\t$tmpPosQ1\t$tmpPosQ3\t$tmpPosIQR\t$tmpPosMAD\n";
    print OUT "\tneg\t$tmpNegMin\t$tmpNegMax\t$tmpNegMedian\t$tmpNegMean\t$tmpNegStdev\t$tmpNegQ1\t$tmpNegQ3\t$tmpNegIQR\t$tmpNegMAD\n";
    print OUT "\n";
    print OUT "\tposBound\t$tmpPosBound\n";
    print OUT "\tnegBound\t$tmpNegBound\n";

    close(OUT);
    
    return($tmpPosBound,$tmpNegBound);
}

sub filterMatrix($$$$) {
    my $matrixObject=shift;
    my $matrix=shift;
    my $anchorOutliers=shift;
    my $manualOutliers=shift;
    
    my $inc2header=$matrixObject->{ inc2header };
    my $numYHeaders=$matrixObject->{ numYHeaders };
    my $numXHeaders=$matrixObject->{ numXHeaders };
    my $missingValue=$matrixObject->{ missingValue };
    my $verbose=$matrixObject->{ verbose };
    
    my $nAnchorOutliers=0;
    my $nManualOutliers=0;
    
    my $nNAs=0;
    
    my %filteredMatrix=();
    for(my $y=0;$y<$numYHeaders;$y++) {
        my $yHeader=$inc2header->{ y }->{$y};
        
        for(my $x=0;$x<$numXHeaders;$x++) {
            my $xHeader=$inc2header->{ x }->{$x};    
            
            my $score=$missingValue;
            $score = $matrix->{$y}->{$x} if(defined($matrix->{$y}->{$x}));
            
            $nAnchorOutliers++ if( (exists($anchorOutliers->{$yHeader})) or (exists($anchorOutliers->{$xHeader})) );
            $nManualOutliers++ if( (exists($manualOutliers->{$yHeader})) or (exists($manualOutliers->{$xHeader})) );
            
            $score = "NA" if( (exists($anchorOutliers->{$yHeader})) or (exists($anchorOutliers->{$xHeader})) );
            $score = "NA" if( (exists($manualOutliers->{$yHeader})) or (exists($manualOutliers->{$xHeader})) );
        
            $nNAs++ if($score eq "NA");
            
            next if( ($score eq $missingValue) or (($score ne "NA") and ($missingValue ne "NA") and ($score == $missingValue)) );
            
            $filteredMatrix{$y}{$x}=$score;
        }    
    }
    
    print STDERR "\tanchorOutliers\t$nAnchorOutliers\n" if($verbose);
    print STDERR "\tmanualOutliers\t$nManualOutliers\n" if($verbose);
    print STDERR "\t$nNAs NA data points\n" if($verbose);
    
    
    return(\%filteredMatrix);
}

my %options;
my $results = GetOptions( \%options,'inputMatrix|i=s','verbose|v','output|o=s','loessObjectFile|lof=s','manualOutlierFile|mof=s','minDistance|minDist=i','maxDistance|maxDist=i','cisAlpha|ca=f','disableIQRFilter|dif','cisApproximateFactor|caf=i','factorMode|fm=s','includeCis|ic','includeTrans|it','iqrMultiplier|im=f','logTransform|lt=f','excludeZero|ez') or croak help();
my ($ret,$inputMatrix,$verbose,$output,$loessObjectFile,$manualOutlierFile,$minDistance,$maxDistance,$cisAlpha,$disableIQRFilter,$cisApproximateFactor,$factorMode,$includeCis,$includeTrans,$iqrMultiplier,$logTransform,$excludeZero) = check_options( \%options );
    
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
my ($matrix)=getData($inputMatrix,$matrixObject,$verbose,$minDistance,$maxDistance);

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
    print STDERR "\tdone\n" if($verbose);
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

print STDERR "generating normalization zScores for Y axis Fragments...\n" if($verbose);
my $yFactorFile=$output.".y.zScores";
my $normalPrimerData={};
$normalPrimerData=getRowColFactor($matrixObject,$inputMatrix,$output,$includeCis,$minDistance,$maxDistance,$includeTrans,$logTransform,$loess,$factorMode,$yFactorFile,$inputMatrixName,$cisApproximateFactor,$excludeZero) if(($includeCis) or ($includeTrans));
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

my $transposedPrimerData={};
if($symmetrical != 1) {
    print STDERR "rotating matrix...\n" if($verbose);
    my $transposedMatrix=transposeMatrix($inputMatrix);
    my $transposedMatrixObject=getMatrixObject($transposedMatrix);
    print STDERR "\tdone\n" if($verbose);

    print STDERR "\n" if($verbose);

    print STDERR "generating normalization zScores for X axis Fragments...\n" if($verbose);
    my $xFactorFile=$output.".x.zScores";
    $transposedPrimerData=getRowColFactor($transposedMatrixObject,$transposedMatrix,$output,$includeCis,$minDistance,$maxDistance,$includeTrans,$logTransform,$loess,$factorMode,$xFactorFile,$inputMatrixName,$cisApproximateFactor,$excludeZero) if(($includeCis) or ($includeTrans));
    print STDERR "\tdone\n" if($verbose);
    
    system("rm ".$transposedMatrix);
    
    print STDERR "\n" if($verbose);
    
}

print STDERR "validating headers...\n" if($verbose);
# validate primer zScore hashes
foreach my $yPrimer (keys %$normalPrimerData) {
    my $yPrimerFactor=$normalPrimerData->{$yPrimer}->{ factor };
    if(exists($transposedPrimerData->{$yPrimer})) {
        next if($yPrimerFactor eq "NA");
        next if($normalPrimerData->{$yPrimer}->{ factor } eq "NA");
        croak "$inputMatrixName - matrix headers symmetrical, but data is not" if($transposedPrimerData->{$yPrimer}->{ factor } != $yPrimerFactor);
    }
}
foreach my $xPrimer (keys %$transposedPrimerData) {
    my $xPrimerFactor=$transposedPrimerData->{$xPrimer}->{ factor };
    if(exists($normalPrimerData->{$xPrimer})) {
        next if($xPrimerFactor eq "NA");
        next if($normalPrimerData->{$xPrimer}->{ factor } eq "NA");
        croak "$inputMatrixName - matrix headers symmetrical, but data is not ($xPrimer)" if($normalPrimerData->{$xPrimer}->{ factor } != $xPrimerFactor);
    }
}
print STDERR "\tvalidated\n" if($verbose);

print STDERR "\n" if($verbose);

# now combine x and y hashes - symmetrical is AOK
my %primerDataHash=();
%primerDataHash=(%$normalPrimerData,%$transposedPrimerData);
my $rowcolData=\%primerDataHash;

# calculate bounds of pos and neg values (cannot assume they will be normal around 0)
# re-scale numbers to be symmetrical -1 -> 0 and 0 -> 1
print STDERR "calculating bounds...\n" if($verbose);
my $tmpPosBound="NA";
my $tmpNegBound="NA";
($tmpPosBound,$tmpNegBound)=calculateBounds($matrixObject,$rowcolData,$iqrMultiplier,$output,$scriptPath,$commentLine) if(keys %{$rowcolData} > 0);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

# write primer factors
print STDERR "writing primer factors\n" if($verbose);
writePrimerFactors($rowcolData,$output,$scriptPath,$commentLine) if(keys %{$rowcolData} > 0);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

print STDERR "writing outliers...\n" if($verbose);
my $anchorOutliers=writeOutliers($matrixObject,$rowcolData,$output,$tmpNegBound,$tmpPosBound,$commentLine) if(keys %{$rowcolData} > 0);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

my $manualOutliers={};
if($manualOutlierFile ne "") {
    print STDERR "reading in manual outlier file...\n" if($verbose);
    $manualOutliers=processManualOutliers($manualOutlierFile);
    print STDERR "\tdone\n" if($verbose);
    print STDERR "\n" if($verbose);
}

print STDERR "filtering matrix...\n" if($verbose);
my $filteredMatrix=filterMatrix($matrixObject,$matrix,$anchorOutliers,$manualOutliers);

print STDERR "\n" if($verbose);

my $filteredMatrixFile=$output.".anchorFiltered.matrix.gz";
print STDERR "writing matrix to file ($filteredMatrixFile)...\n" if($verbose);
writeMatrix($filteredMatrix,$inc2header,$filteredMatrixFile,$missingValue,$commentLine);
print STDERR "\tcomplete\n" if($verbose);

print STDERR "\n" if($verbose);
