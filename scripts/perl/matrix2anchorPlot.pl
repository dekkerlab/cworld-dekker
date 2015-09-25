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

    my ($inputMatrix,$verbose,$output,$loessObjectFile,$cisAlpha,$disableIQRFilter,$minDistance,$maxDistance,$cisApproximateFactor,$excludeZero,$highlightMatrix,$yMin,$yMax,$subsetListFile);
    
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

    if( exists($opts->{ highlightMatrix }) ) {
        $highlightMatrix = $opts->{ highlightMatrix };
    } else {
        $highlightMatrix = "NA";
    }    
    
    if( exists($opts->{ yMin }) ) {
        $yMin = $opts->{ yMin };
    } else {
        $yMin = "NA";
    }    
    
    if( exists($opts->{ yMax }) ) {
        $yMax = $opts->{ yMax };
    } else {
        $yMax = "NA";
    }    
    
    if( exists($opts->{ subsetListFile }) ) {
        $subsetListFile = $opts->{ subsetListFile };
    } else {
        $subsetListFile = "";
    }
    
    $ret->{ inputMatrix }=$inputMatrix;
    $ret->{ verbose }=$verbose;
    $ret->{ output }=$output;
    $ret->{ loessObjectFile }=$loessObjectFile;
    $ret->{ cisAlpha }=$cisAlpha;
    $ret->{ disableIQRFilter }=$disableIQRFilter;
    $ret->{ minDistance }=$minDistance;
    $ret->{ maxDistance }=$maxDistance;
    $ret->{ cisApproximateFactor }=$cisApproximateFactor;
    $ret->{ excludeZero }=$excludeZero;
    $ret->{ highlightMatrix }=$highlightMatrix;
    $ret->{ yMin }=$yMin;
    $ret->{ yMax }=$yMax;
    $ret->{ subsetListFile }=$subsetListFile;
    
    return($ret,$inputMatrix,$verbose,$output,$loessObjectFile,$cisAlpha,$disableIQRFilter,$minDistance,$maxDistance,$cisApproximateFactor,$excludeZero,$highlightMatrix,$yMin,$yMax,$subsetListFile);
}

sub intro() {
    print STDERR "\n";
    
    print STDERR "Tool:\t\t".$tool."\n";
    print STDERR "Version:\t".$cworld::dekker::VERSION."\n";
    print STDERR "Summary:\ttransform each row/col into 4C style 'anchor' plot.\n";
    
    print STDERR "\n";
}

sub help() {
    intro();
    
    print STDERR "Usage: perl matrix2anchorPlot.pl [OPTIONS] -i <inputMatrix>\n";
    
    print STDERR "\n";
    
    print STDERR "Required:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-i", "[]", "input matrix file");
    
    print STDERR "\n";
    
    print STDERR "Options:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-v", "[]", "FLAG, verbose mode");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-o", "[]", "prefix for output file(s)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--lof", "[]", "optional, pre-calculated *.loess file");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ca", "[0.01]", "lowess alpha value, fraction of datapoints to smooth over");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--dif", "[]", "FLAG, disable loess IQR (outlier) filter");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--caf", "[1000]", "cis approximate factor to speed up loess, genomic distance / -caffs");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ez", "[]", "FLAG, ignore 0s in all calculations");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--minDist", "[]", "minimum allowed interaction distance, exclude < N distance (in BP)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--maxDist", "[]", "maximum allowed interaction distance, exclude > N distance (in BP)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--hm", "[]", "highlight matrix file - binary flag for specific interactions");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ymin", "[]", "optional, ymin value for anchor plot y-axis");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ymax", "[]", "optional, ymax value for anchor plot y-axis");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--slf", "[]", "optional, list of headers to draw, all other skipped");
    
    print STDERR "\n";
    
    print STDERR "Notes:";
    print STDERR "
    This script can transform each row/col into 4C style 'anchor' plot.
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

sub processSubsetListFile($) {
    my $subsetListFile=shift;
    
    my %subsetList=();
    
    open(IN,inputWrapper($subsetListFile)) or croak "Could not open file [$subsetListFile] - $!";
    while(my $line = <IN>) {
        chomp($line);
        next if(($line eq "") or ($line =~ m/^#/));
        
        $subsetList{$line}=1;
    }
    close(IN);
    
    return(\%subsetList);
}

sub drawAnchorPlots($$$$$$$$$$$$;$$$) {
    #required
    my $inputMatrix=shift;
    my $output=shift;
    my $includeCis=shift;
    my $minDistance=shift;
    my $maxDistance=shift;
    my $includeTrans=shift;
    my $loess=shift;
    my $anchorPlotFolder=shift;
    my $anchorPlotScriptPath=shift;
    my $subsetList=shift;
    my $yMin=shift;
    my $yMax=shift;
    #optional
    my $cisApproximateFactor=1;
    $cisApproximateFactor=shift if @_;
    my $highlightMatrixFile="NA";
    $highlightMatrixFile=shift if @_;
    my $verbose=0;
    $verbose = shift if @_;
    
    my $cwd = getcwd();
    
    my $matrixObject=getMatrixObject($inputMatrix);
    my $inc2header=$matrixObject->{ inc2header };
    my $header2inc=$matrixObject->{ header2inc };
    my $numYHeaders=$matrixObject->{ numYHeaders };
    my $numXHeaders=$matrixObject->{ numXHeaders };
    
    # get matrix data
    my $matrix={};
    ($matrix)=getData($inputMatrix,$matrixObject);
    
    # process highlight matrix (if given)
    my $highlightMatrix={};
    my $highlightMatrixObject={};
    if($highlightMatrixFile ne "NA") {
        $highlightMatrixObject=getMatrixObject($highlightMatrixFile);
        ($highlightMatrix)=getData($highlightMatrixFile,$highlightMatrixObject);
    }
    
    my $pcComplete=0;    
    for(my $y=0;$y<$numYHeaders;$y++) {
        my $yHeader=$inc2header->{ y }->{$y};
        my $yHeaderObject=getHeaderObject($yHeader);
        
        # if subset list option provided, only process those primers in hash
        next if((keys %{$subsetList} > 0) and (!exists($subsetList->{$yHeader})));
        
        my $yHeaderSubName=$yHeaderObject->{ subName };
        my $anchorSubName=$output."___".$yHeaderSubName;
        my $tmpAnchorPlotFile=$anchorPlotFolder.$anchorSubName.".3C";
        
        if($pcComplete <= 100) {
            print STDERR "\e[A"  if(($verbose) and ($y != 0));
            printf STDERR "\tanchorPlot\t%.2f%% complete ($y/$numYHeaders)...\n", $pcComplete if($verbose);
            $pcComplete = round((($y/$numYHeaders)*100),2);
        }
        
        # clean up any existing files
        system("rm '".$tmpAnchorPlotFile."'") if(-e($tmpAnchorPlotFile));
        system("rm '".$tmpAnchorPlotFile.".expected'") if(-e($tmpAnchorPlotFile.".expected"));
        system("rm '".$tmpAnchorPlotFile.".png'") if(-e($tmpAnchorPlotFile.".png"));
        
        open(OUT,outputWrapper($tmpAnchorPlotFile,"",1)) or croak "Could not open file [$tmpAnchorPlotFile] - $!";
        print OUT "yHeader\txHeader\tinteractionMidpoint\tinteractionDistance\tloessExpected\tloessStdev\tcScore\thighlight\n";
        
        my $drawFlag=0;
        for(my $x=0;$x<$numXHeaders;$x++) {
            my $xHeader=$inc2header->{ x }->{$x};    
            my $xHeaderObject=getHeaderObject($xHeader);
            
            my $xHeaderStart=$xHeaderObject->{ start };
            my $xHeaderEnd=$xHeaderObject->{ end };
            my $xHeaderMidpoint = round(($xHeaderStart + $xHeaderEnd) / 2);
    
            my $interactionDistance=getInteractionDistance($matrixObject,$yHeaderObject,$xHeaderObject,$cisApproximateFactor);
            my $realInteractionDistance=getInteractionDistance($matrixObject,$yHeaderObject,$xHeaderObject,1);
            
            # force classify interaction options - cis only - no limit
            my $interactionClassification=classifyInteraction($matrixObject,$includeCis,$minDistance,$maxDistance,$includeTrans,$yHeaderObject,$xHeaderObject);
                
            my $cScore=$matrixObject->{ missingValue };
            $cScore = $matrix->{$y}->{$x} if(defined($matrix->{$y}->{$x}));
            
            $cScore = "NA" if(($cScore eq "") or ($cScore =~ /^NULL$/i) or ($cScore =~ /^NA$/i) or ($cScore =~ /inf$/i) or ($cScore eq "NA"));
            $cScore = "NA" if($interactionClassification ne "USABLE");
                        
            next if($cScore eq "NA");
            
            #allow draw
            $drawFlag = 1;
            
            croak  "distance ($interactionDistance) does not exist!" if(!exists($loess->{$interactionDistance}));
            
            my $loessExpected = "NA";
            my $loessStdev = "NA";
            $loessExpected = $loess->{$interactionDistance}->{ loess } if(exists($loess->{$interactionDistance}->{ loess }));
            $loessStdev = $loess->{$interactionDistance}->{ stdev } if(exists($loess->{$interactionDistance}->{ stdev }));
            
            my $highlight="NA";
            $highlight=$highlightMatrixObject->{ missingValue } if(defined($highlightMatrixObject->{ missingValue }));
            
            $highlight = $highlightMatrix->{$y}->{$x} if(exists($highlightMatrix->{$y}->{$x}));
            $highlight = "NA" if(($highlight eq "") or ($highlight =~ /^NULL$/i) or ($highlight =~ /^NA$/i) or ($highlight =~ /inf$/i) or ($highlight eq "NA"));
            $highlight = "NA" if($interactionClassification ne "USABLE");
                        
            print OUT "$yHeader\t$xHeader\t$xHeaderMidpoint\t$realInteractionDistance\t$loessExpected\t$loessStdev\t$cScore\t$highlight\n";
            
        }
        
        system("rm '$tmpAnchorPlotFile'") if($drawFlag == 0);
        next if($drawFlag == 0);
        
        my $anchorObject=getHeaderObject($yHeader);
        
        my $anchorShortName=$anchorObject->{ subName };
        my $anchorStart = $anchorObject->{ start };
        my $anchorEnd = $anchorObject->{ end };
        my $anchorMidpoint = round(($anchorStart + $anchorEnd) / 2);
        
        my $expectedPlotFile=$tmpAnchorPlotFile.".expected";
        open(OUT,outputWrapper($expectedPlotFile)) or croak "Could not open file [$expectedPlotFile] - $!";
        print OUT "interactionMidpoint\tloessExpected\tloessStdev\n";
        
        for my $interactionDistance ( sort {$a<=>$b} keys %{$loess}) {
            
            my $loessExpected = "NA";
            my $loessStdev = "NA";
            $loessExpected = $loess->{$interactionDistance}->{ loess } if(exists($loess->{$interactionDistance}->{ loess }));
            $loessStdev = $loess->{$interactionDistance}->{ stdev } if(exists($loess->{$interactionDistance}->{ stdev }));
            
            my $loessUpstreamMidpoint = ($anchorMidpoint - ($interactionDistance*$cisApproximateFactor));
            my $loessDownstreamMidpoint = ($anchorMidpoint + ($interactionDistance*$cisApproximateFactor));
            
            next if(($loessExpected eq "NA") or ($loessStdev eq "NA"));
            
            print OUT "$loessUpstreamMidpoint\t$loessExpected\t$loessStdev\n";
            print OUT "$loessDownstreamMidpoint\t$loessExpected\t$loessStdev\n";
        }
        
        close(OUT);
        
        system("Rscript '".$anchorPlotScriptPath."' '".$cwd."' '".$tmpAnchorPlotFile."' '".$expectedPlotFile."' '".$anchorShortName."' ".$anchorStart." ".$anchorEnd." ".$yMin." ".$yMax." &> /dev/null");
        
        croak "error drawing plot [$tmpAnchorPlotFile]" if(!(-e $tmpAnchorPlotFile.".png"));
        
        system("rm '$tmpAnchorPlotFile'");
        system("rm '$expectedPlotFile'");
    
    }
    
    print STDERR "\e[A" if($verbose);
    printf STDERR "\tanchorPlot\t%.2f%% complete ($numYHeaders/$numYHeaders)...\n", 100 if($verbose);
    
    close(OUT);
    
}

my %options;
my $results = GetOptions( \%options,'inputMatrix|i=s','verbose|v','output|o=s','loessObjectFile|lof=s','cisAlpha|ca=f','disableIQRFilter|dif=s','minDistance|minDist=i','maxDistance|maxDist=i','cisApproximateFactor|caf=i','excludeZero|ez','highlightMatrix|hm=s','yMin|ymin=i','yMax|ymax=i','subsetListFile|slf=s') or croak help();
my ($ret,$inputMatrix,$verbose,$output,$loessObjectFile,$cisAlpha,$disableIQRFilter,$minDistance,$maxDistance,$cisApproximateFactor,$excludeZero,$highlightMatrix,$yMin,$yMax,$subsetListFile)=check_options( \%options );

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

my $anchorPlotFolder="anchor__".$inputMatrixName."/";
system("mkdir -p '$anchorPlotFolder'");

my $subsetList={};
$subsetList=processSubsetListFile($subsetListFile) if(($subsetListFile ne "") and (-e($subsetListFile)));

#validating headers
if($highlightMatrix ne "NA") {
    print STDERR "validating matrices...\n";
    validateIdenticalMatrixStructure($inputMatrix,$highlightMatrix);
    print STDERR "\tdone\n\n";
}

my $excludeCis=0;
my $excludeTrans=1;
my $includeCis=flipBool($excludeCis);
my $includeTrans=flipBool($excludeTrans);

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

# y axis is the anchor
my $anchorPlotScriptPath=$scriptPath."/R/anchorPlot.R";

print STDERR "plotting anchor plots...\n" if($verbose);
drawAnchorPlots($inputMatrix,$output,$includeCis,$minDistance,$maxDistance,$includeTrans,$loess,$anchorPlotFolder,$anchorPlotScriptPath,$subsetList,$yMin,$yMax,$cisApproximateFactor,$highlightMatrix,$verbose) if($includeCis);

print STDERR "\n" if($verbose);

if( !(isSymmetrical($inputMatrix)) )  {

    my $transposedMatrix=transposeMatrix($inputMatrix);
    croak "transposedMatrix [$transposedMatrix] does not exist" if(!(-e $transposedMatrix));
    
    my $transposedHighlightMatrix="NA";
    if($highlightMatrix ne "NA") {
        my $transposedHighlightMatrix=transposeMatrix($highlightMatrix);
        croak "transposedHighlightMatrix [$transposedHighlightMatrix] does not exist" if(!(-e $transposedHighlightMatrix));
    }
    
    print STDERR "plotting anchor plots...\n" if($verbose);
    drawAnchorPlots($transposedMatrix,$output,$includeCis,$minDistance,$maxDistance,$includeTrans,$loess,$anchorPlotFolder,$anchorPlotScriptPath,$subsetList,$yMin,$yMax,$cisApproximateFactor,$transposedHighlightMatrix,$verbose) if($includeCis);
    
    system("rm '$transposedMatrix'");
    system("rm '$transposedHighlightMatrix'") if($transposedHighlightMatrix ne "NA");
    
    print STDERR "\n" if($verbose);
}