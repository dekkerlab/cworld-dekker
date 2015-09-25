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

    my ($inputMatrix_1,$inputMatrix_2,$verbose,$output,$plotPerDistance,$excludeCis,$excludeTrans,$minDistance,$maxDistance,$correlateMode,$excludeZero,$excludeDiagonal,$logTransform,$fixedBinPopulationSize,$outlierFraction,$ymin,$ymax,$xmin,$xmax,$tmpDir);
    
    my $ret={};
    
    if( exists($opts->{ inputMatrix_1 }) ) {
        $inputMatrix_1 = $opts->{ inputMatrix_1 };
    } else {
        print STDERR "\nERROR: Option inputMatrix_1|1 is required.\n";
        help();
    }
    
    if( exists($opts->{ inputMatrix_2 }) ) {
        $inputMatrix_2 = $opts->{ inputMatrix_2 };
    } else {
        print STDERR "\nERROR: Option inputMatrix_1|1 is required.\n";
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
    
    if( exists($opts->{ plotPerDistance })) {
        $plotPerDistance=1;
    } else {
        $plotPerDistance=0;
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
    
    if( exists($opts->{ maxDistance }) ) {
        $maxDistance = $opts->{ maxDistance };
    } else {
        $maxDistance = undef;
    }
    
    if( exists($opts->{ minDistance }) ) {
        $minDistance = $opts->{ minDistance };
    } else {
        $minDistance=undef;
    }
    
    if( exists($opts->{ correlateMode }) ) {
        $correlateMode = $opts->{ correlateMode };
    } else {
        $correlateMode = "pearson";
    }
    
    if( exists($opts->{ excludeZero }) ) {
        $excludeZero = 1;
    } else {
        $excludeZero = 0;
    }
    
    if( exists($opts->{ excludeDiagonal }) ) {
        $excludeDiagonal = 1;
    } else {
        $excludeDiagonal = 0;
    }
    
    if( exists($opts->{ logTransform }) ) {
        $logTransform = $opts->{ logTransform };
    } else {
        $logTransform = 0;
    }
    
    if( exists($opts->{ fixedBinPopulationSize }) ) {
        $fixedBinPopulationSize = $opts->{ fixedBinPopulationSize };
    } else {
        $fixedBinPopulationSize = 0;
    }    
    
    if( exists($opts->{ outlierFraction }) ) {
        $outlierFraction = $opts->{ outlierFraction };
    } else {
        $outlierFraction = 0.05;
    }
    
    if( exists($opts->{ ymin }) ) {
        $ymin = $opts->{ ymin };
    } else {
        $ymin = "NA";
    }
    
    if( exists($opts->{ ymax }) ) {
        $ymax = $opts->{ ymax };
    } else {
        $ymax = "NA";
    }
    
    if( exists($opts->{ xmin }) ) {
        $xmin = $opts->{ xmin };
    } else {
        $xmin = "NA";
    }
    
    if( exists($opts->{ xmax }) ) {
        $xmax = $opts->{ xmax };
    } else {
        $xmax = "NA";
    }
    
    if( exists($opts->{ tmpDir }) ) {
        $tmpDir = $opts->{ tmpDir };
        $tmpDir =~ s/\/$//;
    } else {
        $tmpDir = "/tmp";
    }
    
    $ret->{ inputMatrix_1 }=$inputMatrix_1;
    $ret->{ inputMatrix_2 }=$inputMatrix_2;
    $ret->{ verbose }=$verbose;
    $ret->{ output }=$output;
    $ret->{ plotPerDistance }=$plotPerDistance;
    $ret->{ excludeCis }=$excludeCis;
    $ret->{ excludeTrans }=$excludeTrans;
    $ret->{ minDistance }=$minDistance;
    $ret->{ maxDistance }=$maxDistance;
    $ret->{ correlateMode }=$correlateMode;
    $ret->{ excludeZero }=$excludeZero;
    $ret->{ excludeDiagonal }=$excludeDiagonal;
    $ret->{ logTransform }=$logTransform;
    $ret->{ fixedBinPopulationSize }=$fixedBinPopulationSize;
    $ret->{ outlierFraction }=$outlierFraction;
    $ret->{ ymin }=$ymin;
    $ret->{ ymax }=$ymax;
    $ret->{ xmin }=$xmin;
    $ret->{ xmax }=$xmax;
    $ret->{ tmpDir }=$tmpDir;
    
    return($ret,$inputMatrix_1,$inputMatrix_2,$verbose,$output,$plotPerDistance,$excludeCis,$excludeTrans,$minDistance,$maxDistance,$correlateMode,$excludeZero,$excludeDiagonal,$logTransform,$fixedBinPopulationSize,$outlierFraction,$ymin,$ymax,$xmin,$xmax,$tmpDir);
}

sub intro() {
    print STDERR "\n";
    
    print STDERR "Tool:\t\t".$tool."\n";
    print STDERR "Version:\t".$cworld::dekker::VERSION."\n";
    print STDERR "Summary:\tperforms correlation between two matrices\n";
    
    print STDERR "\n";
}

sub help() {
    intro();
    
    print STDERR "Usage: perl correlateMatrices.pl [OPTIONS] -i1 <inputMatrix_1> -i2 <inputMatrix_2>\n";
    
    print STDERR "\n";
    
    print STDERR "Required:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-1", "[]", "input matrix 1 file");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-2", "[]", "input matrix 2 file");
    
    print STDERR "\n";
    
    print STDERR "Options:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-v", "[]", "FLAG, verbose mode");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-o", "[]", "prefix for output file(s)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ppd", "[]", "FLAG, plot correlation per distance");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ec", "[]", "FLAG, exclude CIS data");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--et", "[]", "FLAG, exclude TRANS data");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ez", "[]", "FLAG, exclude zeros from all calculations");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ed", "[]", "FLAG, exclude diagonal bin (y=x) from all calculations");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--minDist", "[]", "minimum allowed interaction distance, exclude < N distance (in BP)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--maxDist", "[]", "maximum allowed interaction distance, exclude > N distance (in BP)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--cm", "[]", "optional, correlation mode (pearson,spearman)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--lt", "[]", "optional, log transform data (2,10)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--fbp", "[]", "optional, fixed bin population size (num of interactions) (-1=auto, 0= ff)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--of", "[]", "optional, outlierFraction [0-1], removed from top/bottom during correlation, black=outliers");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ymin", "[]", "optional, y min value of plot");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ymax", "[]", "optional, y max value of plot");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--xmin", "[]", "optional, x min value of plot");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--xmax", "[]", "optional, x max value of plot");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--tmp", "[/tmp]", "optional, tmp direction for tmp files");
    
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

sub correlatePerDistance($$$$) {
    my $pairwiseFile=shift;
    my $scriptPath=shift;
    my $output=shift;
    my $tmpDir=shift;
    
    $tmpDir=createTmpDir($tmpDir);
    
    my $correlateList=$scriptPath."/R/correlateList.R";
    
    my $correlationPerDistanceFile=$output.".correlationPerDistance.txt.gz";
    open(OUT,outputWrapper($correlationPerDistanceFile)) or croak "Could not open file [$correlationPerDistanceFile] - $!";
    print OUT "interactionDistance\tcorrelation\n";
    
    my $tmpFile="";
    my $lastInteractionDistance="NA";
    my $TMP;
    
    open(IN,inputWrapper($pairwiseFile)) or croak "Could not open file [$pairwiseFile] - $!";
    while(my $line = <IN>) {
        chomp($line);
        next if(($line eq "") or ($line =~ m/^#/));
        
        my ($interactionDistance,$score_1,$score_2)=split(/\t/,$line);
        
        if($lastInteractionDistance eq "NA") {
            $tmpFile=$tmpDir."tmp_".$interactionDistance.".txt.gz";
            open($TMP,outputWrapper($tmpFile)) or croak "Could not open file [$tmpFile] - $!";
            print $TMP "interactionDistance\tscore_1\tscore_2\n";
            $lastInteractionDistance=$interactionDistance;
        }
        
        if($interactionDistance != $lastInteractionDistance) {
            close($TMP);
            my $r=`Rscript $correlateList $tmpFile`;
            print OUT "$lastInteractionDistance\t$r\n";
            
            $tmpFile=$tmpDir."tmp_".$interactionDistance.".txt.gz";
            open($TMP,outputWrapper($tmpFile)) or croak "Could not open file [$tmpFile] - $!";
            print $TMP "interactionDistance\tscore_1\tscore_2\n";
            
        } else {
            print $TMP "$interactionDistance\t$score_1\t$score_2\n";
        }
        
        $lastInteractionDistance=$interactionDistance;
    }
    close($TMP);
    my $r=`Rscript $correlateList $tmpFile`;
    print OUT "$lastInteractionDistance\t$r\n";
            
    close(IN);
    
    close(OUT);
    
    system("Rscript ".$scriptPath."correlateMatrices/scripts/correlateMatricesPerDistance.R $correlationPerDistanceFile");
    
    removeTmpDir($tmpDir);
}

sub fixBinPopulation($$$) {
    my $pairwiseFile=shift;
    my $output=shift;
    my $distancePopulation_size=shift;
    
    my $binIndexFile=$output.".sorted.binIndex.txt.gz";
    open(OUT,outputWrapper($binIndexFile)) or croak "Could not open file [$binIndexFile] - $!";
    
    my %binIndex2distance=();
    
    my $binIndex=0;
    my $binPopulation=0;
    my $lastInteractionDistance="NA";
    my $binMin="NA";
    my $binMax="NA";
    
    open(IN,inputWrapper($pairwiseFile)) or croak "Could not open file [$pairwiseFile] - $!";
    while(my $line = <IN>) {
        chomp($line);
        next if(($line eq "") or ($line =~ m/^#/));
        
        my ($interactionDistance,$score_1,$score_2)=split(/\t/,$line);
        
        # skip all distance < 0
        print OUT "$interactionDistance\t$score_1\t$score_2\n" if($interactionDistance < 0);
        next if($interactionDistance < 0);
        
        if($lastInteractionDistance eq "NA") {
            $lastInteractionDistance=$interactionDistance;
            $binMin=$interactionDistance;
        }
        
        if(($interactionDistance != $lastInteractionDistance) and ($binPopulation > $distancePopulation_size)) {
            my $binMidpoint=(($binMin+$binMax)/2);
            print STDERR "$binIndex\t$binPopulation\t$binMin - $binMax\t$binMidpoint\n";
            $binIndex2distance{$binIndex}=$binMidpoint;
            
            $binIndex++;
            $binPopulation=1;
            $binMin=$interactionDistance;
            $binMax=$interactionDistance;
        } else {
            $binPopulation++;
        }
        $lastInteractionDistance=$interactionDistance;
        $binMax=$interactionDistance;
        
        print OUT "$binIndex\t$score_1\t$score_2\n";
    }
    close(IN);
    
    close(OUT);
    my $binMidpoint=(($binMin+$binMax)/2);
    print STDERR "$binIndex\t$binPopulation\t$binMin - $binMax\t$binMidpoint\n";
    $binIndex2distance{$binIndex}=$binMidpoint;
    
    # now replace bin indeces with bin midpoints
    
    my $binMidpointFile=$output.".sorted.binMidpoint.txt.gz";
    open(OUT,outputWrapper($binMidpointFile)) or croak "Could not open file [$binMidpointFile] - $!";
    
    open(IN,inputWrapper($binIndexFile)) or croak "Could not open file [$binIndexFile] - $!";
    while(my $line = <IN>) {
        chomp($line);
        next if(($line eq "") or ($line =~ m/^#/));
        
        my ($interactionDistance,$score_1,$score_2)=split(/\t/,$line);
        
        $interactionDistance=$binIndex2distance{$interactionDistance} if(exists($binIndex2distance{$interactionDistance}));
        
        print OUT "$interactionDistance\t$score_1\t$score_2\n";
        
    }
    close(IN);
    
    close(OUT);
    
    return($binMidpointFile);
    
}
        

my %options;
my $results = GetOptions( \%options,'inputMatrix_1|1=s','inputMatrix_2|2=s','verbose|v','output|o=s','plotPerDistance|ppd','excludeCis|ec','excludeTrans|et','minDistance|minDist=i','maxDistance|maxDist=i','correlateMode|cm=s','excludeZero|ez','excludeDiagonal|ed','logTransform|lt=f','fixedBinPopulationSize|fbps=s','outlierFraction|of=i','ymin=i','ymax=i','xmin=i','xmax=i','tmpDir|tmp=s') or croak help();
my ($ret,$inputMatrix_1,$inputMatrix_2,$verbose,$output,$plotPerDistance,$excludeCis,$excludeTrans,$minDistance,$maxDistance,$correlateMode,$excludeZero,$excludeDiagonal,$logTransform,$fixedBinPopulationSize,$outlierFraction,$ymin,$ymax,$xmin,$xmax,$tmpDir)=check_options( \%options );

intro() if($verbose);

#get the absolute path of the script being used
my $cwd = getcwd();
my $fullScriptPath=abs_path($0);
my @fullScriptPathArr=split(/\//,$fullScriptPath);
@fullScriptPathArr=@fullScriptPathArr[0..@fullScriptPathArr-3];
my $scriptPath=join("/",@fullScriptPathArr);
my $commentLine=getScriptOpts($ret,$tool);

# ensure input matrices exist
croak "inputMatrix [$inputMatrix_1] does not exist" if(!(-e $inputMatrix_1));
croak "inputMatrix [$inputMatrix_2] does not exist" if(!(-e $inputMatrix_2));
croak "inputMatrix_1 == inputMatrix_2" if($inputMatrix_1 eq $inputMatrix_2);

my $inputMatrixName_1=getFileName($inputMatrix_1);
my $inputMatrixName_2=getFileName($inputMatrix_2);
$output = $inputMatrixName_1."___".$inputMatrixName_2 if($output eq "");

#validating headers
print STDERR "validating headers ...\n" if($verbose);
validateIdenticalMatrixStructure($inputMatrix_1,$inputMatrix_2);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

# get matrix information
my $matrixObject_1=getMatrixObject($inputMatrix_1,$output,$verbose);
my $inc2header_1=$matrixObject_1->{ inc2header };
my $header2inc_1=$matrixObject_1->{ header2inc };
my $numYHeaders_1=$matrixObject_1->{ numYHeaders };
my $numXHeaders_1=$matrixObject_1->{ numXHeaders };

my $matrix_1={};
($matrix_1)=getData($inputMatrix_1,$matrixObject_1,$verbose,$minDistance,$maxDistance,$excludeCis,$excludeTrans);

print STDERR "\n" if($verbose); 

# get matrix information
my $matrixObject_2=getMatrixObject($inputMatrix_2,$output,$verbose);
my $inc2header_2=$matrixObject_2->{ inc2header };
my $header2inc_2=$matrixObject_2->{ header2inc };
my $numYHeaders_2=$matrixObject_2->{ numYHeaders };
my $numXHeaders_2=$matrixObject_2->{ numXHeaders };

my $matrix_2={};
($matrix_2)=getData($inputMatrix_2,$matrixObject_2,$verbose,$minDistance,$maxDistance,$excludeCis,$excludeTrans);

print STDERR "\n" if($verbose);

# assume both 1 and 2 are equal
my $inc2header=$inc2header_1;
my $header2inc=$header2inc_1;

# do user specified compare mode
print STDERR "dumping pairwise file ...\n" if($verbose);
my ($pairwiseFile,$distHist)=correlateMatrices($output,$matrixObject_1,$matrixObject_2,$matrix_1,$matrix_2,$inc2header,$excludeZero,$excludeDiagonal,$logTransform);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

print STDERR "calculating populations per distance...\n" if($verbose);
my @distArr=();
foreach my $dist (keys %$distHist) {
    my $population=$distHist->{$dist};
    push(@distArr,$population);
}
my $distArrStats=listStats(\@distArr);
my $distancePopulation_max=$distArrStats->{ max };
my $distancePopulation_stdev=$distArrStats->{ stdev };
my $distancePopulation_mean=$distArrStats->{ mean };
my $distancePopulation_median=$distArrStats->{ median };
my $distancePopulation_iqrMean=$distArrStats->{ iqrMean };
print STDERR "\tmax = $distancePopulation_max\n" if($verbose);
print STDERR "\tstdev = $distancePopulation_stdev\n" if($verbose);
print STDERR "\tmean = $distancePopulation_mean\n" if($verbose);
print STDERR "\tmedian = $distancePopulation_median\n" if($verbose);
print STDERR "\tiqrMean = $distancePopulation_iqrMean\n" if($verbose);

print STDERR "\n" if($verbose);

print STDERR "plotting correlation ...\n" if($verbose);
system("Rscript '".$scriptPath."/R/correlateMatrices.R' '".$cwd."' '".$pairwiseFile."' '".$inputMatrixName_1."' '".$inputMatrixName_2."' '".$correlateMode."' ".$outlierFraction." ".$ymin." ".$ymax." ".$xmin." ".$xmax." &> /dev/null");

print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

print STDERR "sorting pairwise file by distance ...\n" if($verbose);
my $sortedPairwiseFile=$output.".correlate.sorted.txt.gz";
system("gunzip -c $pairwiseFile | tail -n +2 | sort -k1,1n | gzip > $sortedPairwiseFile");
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

if($plotPerDistance) {
    
    if($fixedBinPopulationSize != 0) {
        print STDERR "re-binning distances into fixed sizes bins...\n" if($verbose);
        my $desiredPopulationSize=($distancePopulation_max*3);
        $desiredPopulationSize=$fixedBinPopulationSize if($fixedBinPopulationSize > 0);
        print STDERR "\tdesired bin population = $desiredPopulationSize\n" if($verbose);
        $sortedPairwiseFile=fixBinPopulation($sortedPairwiseFile,$output,$desiredPopulationSize);
        print STDERR "\tdone\n" if($verbose);
        print STDERR "\n" if($verbose);
    }

    print STDERR "correlating by distance\n" if($verbose);
    correlatePerDistance($sortedPairwiseFile,$scriptPath,$output,$tmpDir);
    print STDERR "done\n" if($verbose);

    print STDERR "\n" if($verbose);

}

