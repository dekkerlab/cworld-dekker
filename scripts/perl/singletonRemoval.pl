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
    
    my ($inputMatrix,$verbose,$output,$loessObjectFile,$includeCis,$includeTrans,$manualOutlierFile,$cisAlpha,$disableIQRFilter,$minDistance,$maxDistance,$cisApproximateFactor,$cisTrimAmount,$transTrimAmount,$logTransform,$excludeZero);
    
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

    if( exists($opts->{ includeTrans }) ) {
        $includeTrans=1;
    } else {
        $includeTrans=0;
    }
    
    if( exists($opts->{ manualOutlierFile }) ) {
        $manualOutlierFile = $opts->{ manualOutlierFile };
    } else {
        $manualOutlierFile="";
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
    
    if( exists($opts->{ cisApproximateFactor }) ) {
        $cisApproximateFactor = $opts->{ cisApproximateFactor };
    } else {
        $cisApproximateFactor=1000;
    }
    
    if( exists($opts->{ cisTrimAmount }) ) {
        $cisTrimAmount = $opts->{ cisTrimAmount };
    } else {
        $cisTrimAmount = 12;
    }    
    
    if( exists($opts->{ transTrimAmount }) ) {
        $transTrimAmount = $opts->{ transTrimAmount };
    } else {
        $transTrimAmount = 12;
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
    
    if(($manualOutlierFile eq "") and (($includeCis + $includeTrans) == 0)) {
        print STDERR "\nERROR: must choose either to include CIS (--ic) data and/or TRANS (--it) data.\n";
        help();
    }
    
    return($inputMatrix,$verbose,$output,$loessObjectFile,$includeCis,$includeTrans,$manualOutlierFile,$cisAlpha,$disableIQRFilter,$minDistance,$maxDistance,$cisApproximateFactor,$cisTrimAmount,$transTrimAmount,$logTransform,$excludeZero);
}

sub intro() {
    print STDERR "\n";
    
    print STDERR "Tool:\t\tsingletonRemoval.pl\n";
    print STDERR "Version:\t".$cworld::dekker::VERSION."\n";
    print STDERR "Summary:\tdetect and remove singleton outliers\n";
    
    print STDERR "\n";
}

sub help() {
    intro();
    
    print STDERR "Usage: perl singletonRemoval.pl [OPTIONS] -i <inputMatrix>\n";
    
    print STDERR "\n";
    
    print STDERR "Required:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-i", "[]", "input matrix file");  
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ic", "[]", "model CIS data to detect outlier row/cols");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--it", "[]", "model TRANS data to detect outlier row/cols");
    print STDERR "\n";
    
    print STDERR "Options:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-v", "[]", "FLAG, verbose mode");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-o", "[]", "prefix for output file(s)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--lof", "[]", "optional loess object file (pre-calculated loess)");   
    printf STDERR ("\t%-10s %-10s %-10s\n", "--mof", "[]", "optional manual outlier file, 1 header per line to be filtered.");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ca", "[0.01]", "lowess alpha value, fraction of datapoints to smooth over");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--dif", "[]", "FLAG, disable loess IQR (outlier) filter");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--minDist", "[]", "minimum allowed interaction distance, exclude < N distance (in BP)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--maxDist", "[]", "maximum allowed interaction distance, exclude > N distance (in BP)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--caf", "[1000]", "cis approximate factor to speed up loess, genomic distance / -caffs");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--cta", "[12]", "z-score threshold for cis interactions");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--tta", "[12]", "z-score threshold for trans interactions");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ez", "[]", "FLAG, ignore 0s in all calculations");
    print STDERR "\n";
    
    print STDERR "Notes:";
    print STDERR "
    This script can detect/remove singleton outlier interactions.
    Matrix should have row/col headers in the standard my5C format.
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

sub detectOutliers($$$$$) {
    my $matrixFile=shift;
    my $output=shift;
    my $cisTrimAmount=shift;
    my $transTrimAmount=shift;
    my $verbose=shift;
    
    my $cwd = getcwd();
        
    # get matrix information
    my $matrixObject=getMatrixObject($matrixFile,$output);
    my $inc2header=$matrixObject->{ inc2header };
    my $header2inc=$matrixObject->{ header2inc };
    my $numYHeaders=$matrixObject->{ numYHeaders };
    my $numXHeaders=$matrixObject->{ numXHeaders };
    my $missingValue=$matrixObject->{ missingValue };
    my $symmetrical=$matrixObject->{ symmetrical };
    my $equalHeaderFlag=$matrixObject->{ equalHeaderFlag };
    my $inputMatrixName=$matrixObject->{ inputMatrixName };
    $output=$matrixObject->{ output };

    #read Matrix
    my $matrix={};
    ($matrix)=getData($matrixFile,$matrixObject,$verbose);
    
    my %outliers=();
    
    for(my $y=0;$y<$numYHeaders;$y++) {
        my $yHeader=$inc2header->{ y }->{$y};
        my $yHeaderObject=getHeaderObject($yHeader);
        for(my $x=0;$x<$numXHeaders;$x++) {
            my $xHeader=$inc2header->{ x }->{$x};
            my $xHeaderObject=getHeaderObject($xHeader);
        
            my $cScore=$matrixObject->{ missingValue };
            $cScore=$matrix->{$y}->{$x} if(defined($matrix->{$y}->{$x}));
            next if($cScore eq "NA");
            
            my $interactionDistance=getInteractionDistance($matrixObject,$yHeaderObject,$xHeaderObject,1);
            
            next if(($interactionDistance == -1) and ($cScore < $transTrimAmount) and ($cScore > -$transTrimAmount));
            next if(($interactionDistance != -1) and ($cScore < $cisTrimAmount) and ($cScore > -$cisTrimAmount));
            
            # outlier
            
            my $key = $yHeader."___".$xHeader;
            my $type="cis";
            $type="trans" if($interactionDistance == -1);
            
            $outliers{$key}{ score }=$cScore;
            $outliers{$key}{ type }=$type;
        }
    }
    
    return(\%outliers);
        
}

sub filterMatrix($$$$) {
    my $matrixObject=shift;
    my $matrix=shift;
    my $outliers=shift;
    my $manualOutliers=shift;
    
    my $inc2header=$matrixObject->{ inc2header };
    my $header2inc=$matrixObject->{ header2inc };
    my $numYHeaders=$matrixObject->{ numYHeaders };
    my $numXHeaders=$matrixObject->{ numXHeaders };
    my $symmetrical=$matrixObject->{ symmetrical };
    my $verbose=$matrixObject->{ verbose };
    
    my $nOutliers=0;
    my $nManualOutliers=0;
    my $nNAs=0;
    
    my %filteredMatrix=();
    for(my $y=0;$y<$numYHeaders;$y++) {
        my $yHeader=$inc2header->{ y }->{$y};
        
        for(my $x=0;$x<$numXHeaders;$x++) {
            my $xHeader=$inc2header->{ x }->{$x};    
            
            my $score=$matrixObject->{ missingValue };
            $score = $matrix->{$y}->{$x} if(defined($matrix->{$y}->{$x}));
            
            my $interactionKey_above=$yHeader."___".$xHeader;
            my $interactionKey_below=$xHeader."___".$yHeader;
            
            $nOutliers++ if(exists($outliers->{$interactionKey_above}));            
            $nManualOutliers++ if(exists($manualOutliers->{$interactionKey_above}));
            
            $score = "NA" if(exists($outliers->{$interactionKey_above}));
            $score = "NA" if(exists($manualOutliers->{$interactionKey_above}));
            
            # only if symmetrical 
            if($symmetrical) {
                $nOutliers++ if(exists($outliers->{$interactionKey_below}));
                $nManualOutliers++ if(exists($manualOutliers->{$interactionKey_below}));
                
                $score = "NA" if(exists($outliers->{$interactionKey_below}));
                $score = "NA" if(exists($manualOutliers->{$interactionKey_below}));
            }
            
            $nNAs++ if($score eq "NA");
            
            $filteredMatrix{$y}{$x}=$score;
        }    
    }
    
    print STDERR "\toutliers\t$nOutliers\n" if($verbose);
    print STDERR "\tmanualOutliers\t$nManualOutliers\n" if($verbose);
    print STDERR "\t$nNAs NA data points\n" if($verbose);
    
    return(\%filteredMatrix);
}

sub writeOutliers($$$$) {
    my $toRemoveFile=shift;
    my $inputMatrixName=shift;
    my $outliers=shift;
    my $manualOutliers=shift;
    
    
    # output all outliers

    my $time = getDate();

    open(OUT,outputWrapper($toRemoveFile)) or croak "Could not open file [$toRemoveFile] - $!";

    # print the cis outliers
    foreach my $key ( keys %$outliers ) {
        my $score=$outliers->{$key}->{ score };
        my $type=$outliers->{$key}->{ type };
        print OUT "$type\t$key\t$score\n";
    }
    close(OUT);
    
}

my %options;
my $results = GetOptions( \%options,'inputMatrix|i=s','verbose|v','output|o=s','loessObjectFile|lof=s','includeCis|ic','includeTrans|it','manualOutlierFile|mof=i','cisAlpha|ca=f','disableIQRFilter|dif','minDistance|minDist=i','maxDistance|maxDist=i','cisApproximateFactor|caf=i','cisTrimAmount|cta=f','transTrimAmount|tta=f','logTransform|lt=f','excludeZero|ez') or croak help();

my ($inputMatrix,$verbose,$output,$loessObjectFile,$includeCis,$includeTrans,$manualOutlierFile,$cisAlpha,$disableIQRFilter,$minDistance,$maxDistance,$cisApproximateFactor,$cisTrimAmount,$transTrimAmount,$logTransform,$excludeZero)=check_options( \%options );

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

#get zScore matrix
print STDERR "calculating z-score matrix...\n" if($verbose);
my $zScoreMatrix=calculateZscore($matrixObject,$matrix,$loess,$cisApproximateFactor,$excludeZero);
print STDERR "\tdone\n" if($verbose);
 
print STDERR "\n" if($verbose);

 my $zScoreFile=$output.".zScore.matrix.gz";
print STDERR "writing matrix to file ($zScoreFile)...\n" if($verbose);
writeMatrix($zScoreMatrix,$inc2header,$zScoreFile,"NA");
print STDERR "\tcomplete\n" if($verbose);

print STDERR "\n" if($verbose);

print STDERR "detecting outliers ...\n" if($verbose);
my $outliers=detectOutliers($zScoreFile,$output,$cisTrimAmount,$transTrimAmount,$verbose);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

my $manualOutliers={};
if($manualOutlierFile ne "") {
    print STDERR "reading in manual outlier file...\n" if($verbose);
    $manualOutliers=processManualOutliers($manualOutlierFile);
    print STDERR "\tdone\n" if($verbose);
    print STDERR "\n" if($verbose);
}

my $toRemoveFile=$output.".toRemove.txt";
print STDERR "writing $toRemoveFile...\n" if($verbose);
writeOutliers($toRemoveFile,$inputMatrixName,$outliers,$manualOutliers);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

print STDERR "filtering matrix...\n" if($verbose);
my $filteredMatrix=filterMatrix($matrixObject,$matrix,$outliers,$manualOutliers);

print STDERR "\n" if($verbose);

my $filteredMatrixFile=$output.".outlierFiltered.matrix.gz";
print STDERR "writing matrix to file ($filteredMatrixFile)...\n" if($verbose);
writeMatrix($filteredMatrix,$inc2header,$filteredMatrixFile);
print STDERR "\tcomplete\n" if($verbose);

print STDERR "\n" if($verbose);