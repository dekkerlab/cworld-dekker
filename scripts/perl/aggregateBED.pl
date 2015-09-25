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

    my ($inputBedFile,$verbose,$output,$contigBoundFile,$assembly,$windowSize,$windowStep,$windowMode,$missingScore,$cisMode,$excludeZero,$yBound);
    
    my $ret={};
    
    if( exists($opts->{ inputBedFile }) ) {
        $inputBedFile = $opts->{ inputBedFile };
    } else {
        print STDERR "\nERROR: Option inputBedFile|i is required.\n";
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
    
    if( exists($opts->{ contigBoundFile }) ) {
        $contigBoundFile = $opts->{ contigBoundFile };
    } else {
        $contigBoundFile="";
    }
    
    if( exists($opts->{ assembly }) ) {
        $assembly = $opts->{ assembly };
    } else {
        print STDERR "\nERROR: Option assembly|a is required.\n";
        help();
    }
    
    if( exists($opts->{ windowSize }) ) {
        $windowSize = $opts->{ windowSize };
    } else {
        $windowSize = 100000;
    }
    
    if( exists($opts->{ windowStep }) ) {
        $windowStep = $opts->{ windowStep };
    } else {        
        $windowStep=1;
    }    
    
    if( exists($opts->{ windowMode }) ) {
        $windowMode = $opts->{ windowMode };
    } else {
        $windowMode="sum";
    }
    
    if( exists($opts->{ missingScore }) ) {
        $missingScore = $opts->{ missingScore };
    } else {
        $missingScore="NA";
    }

    if( exists($opts->{ cisMode }) ) {
        $cisMode = 1;
    } else {
        $cisMode = 0;
    }
    
    if( exists($opts->{ excludeZero }) ) {
        $excludeZero = 1;
    } else {
        $excludeZero = 0;
    }
    
    if( exists($opts->{ yBound }) ) {
        $yBound = $opts->{ yBound };
    } else {
        $yBound = "";
    }
    
    $ret->{ inputBedFile }=$inputBedFile;
    $ret->{ verbose }=$verbose;
    $ret->{ output }=$output;
    $ret->{ contigBoundFile }=$contigBoundFile;
    $ret->{ assembly }=$assembly;
    $ret->{ windowSize }=$windowSize;
    $ret->{ windowStep }=$windowStep;
    $ret->{ windowMode }=$windowMode;
    $ret->{ missingScore }=$missingScore;
    $ret->{ cisMode }=$cisMode;
    $ret->{ excludeZero }=$excludeZero;
    $ret->{ yBound }=$yBound;

    return($ret,$inputBedFile,$verbose,$output,$contigBoundFile,$assembly,$windowSize,$windowStep,$windowMode,$missingScore,$cisMode,$excludeZero,$yBound);
}

sub intro() {
    print STDERR "\n";
    
    print STDERR "Tool:\t\t".$tool."\n";
    print STDERR "Version:\t".$cworld::dekker::VERSION."\n";
    print STDERR "Summary:\tsliding window aggregrate of BED5 data\n";
    
    print STDERR "\n";
}

sub help() {
    intro();
    
    print STDERR "Usage: perl aggregrateBED.pl [OPTIONS] -i <inputBEDFile> -a <assembly>\n";
    
    print STDERR "\n";
        
    print STDERR "Required:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-i", "[]", "input BED5 file");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--a", "[]", "assembly of BED5 file");
    
    print STDERR "\n";

    print STDERR "Options:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-v", "[]", "FLAG, verbose mode");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-o", "[]", "prefix for output file(s)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--cbf", "[]", "contig bound file (bed3+), to specify start/end of regions/contigs.");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--wsize", "[100000]", "window size, in bp");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--wstep", "[1]", "window step factor. 1 = non overlapping windows, 2-inf = overlapping windows");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--wmode", "[sum]", "window signal aggregrate mode.  available modes are sum,mean,median,binary,min,max,stdev,variance etc");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ms", "[NA]", "suplment missing scores with supplied value. (NA,0)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--cis", "[]", "enable cis mode - seperate track for each contig");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ez", "[]", "FLAG, ignore 0s in all calculations");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--yl", "[-1]", "optional ylimit for plot, by default, plots are autoScaled");
    
    print STDERR "\n";
    
    print STDERR "Notes:";
    print STDERR "
    This program can window a BED5 formatted file into fixed size intervals.
    BED signal is then aggregrated (mean,median,sum etc)\n";
    
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

sub loadContigBoundFile($$$$) {
    my $contigBoundFile=shift;
    my $windowSize=shift;
    my $windowStep=shift;
    my $verbose=shift;
    
    print STDERR "loading contig bounds ...\n" if($verbose);
    
    validateBED($contigBoundFile);
    
    my %contigInfo=();
    my %windowSizes=();
    
    open(IN,inputWrapper($contigBoundFile)) or croak "Could not open file [$contigBoundFile] - $!";
    while(my $line = <IN>) {
        chomp($line);
        next if(($line eq "") or ($line =~ m/^#/));
        
        # skip possible BED headers
        next if($line =~ /^#/);
        next if($line =~ /^track/);
        next if($line =~ /^chrom/);
        
        my @tmp=split(/\t/,$line);
        
        my $chromosome=$tmp[0];
        my $start=$tmp[1];
        my $end=$tmp[2];
        
        $contigInfo{$chromosome}{ start }=$start;
        $contigInfo{$chromosome}{ end }=$end;
            
        my $minIndex=floor($start/($windowSize/$windowStep));
        my $maxIndex=floor($end/($windowSize/$windowStep))+1;
        $windowSizes{$chromosome}{ min } = 0;
        $windowSizes{$chromosome}{ max } = ($maxIndex-$minIndex);
        
        print STDERR "\t$chromosome\t$start\t$end\n" if($verbose);
        
    }
    close(IN);
    
    return(\%contigInfo,\%windowSizes);
    
    print STDERR "\n" if($verbose);
    
}

my %options;
my $results = GetOptions( \%options,'inputBedFile|i=s','verbose|v','output|o=s','contigBoundFile|cbf=s','assembly|a=s','windowSize|wsize=i','windowStep|wstep=i','windowMode|wmode=s','missingScore|ms=s','cisMode|cis','excludeZero|ez','yBound|yb=f') or croak help();
my ($ret,$inputBedFile,$verbose,$output,$contigBoundFile,$assembly,$windowSize,$windowStep,$windowMode,$missingScore,$cisMode,$excludeZero,$yBound)=check_options( \%options );

intro() if($verbose);

#get the absolute path of the script being used
my $cwd = getcwd();
my $fullScriptPath=abs_path($0);
my @fullScriptPathArr=split(/\//,$fullScriptPath);
@fullScriptPathArr=@fullScriptPathArr[0..@fullScriptPathArr-3];
my $scriptPath=join("/",@fullScriptPathArr);

my $windowStepSize=ceil($windowSize/$windowStep);

$output=removeFileExtension(getFileName($inputBedFile))."__smoothed" if($output eq "");

croak "inputFile ($inputBedFile) does not exist!" if(!(-e $inputBedFile));

validateBED($inputBedFile);

my $contigInfo={};
my $windowSizes={};
($contigInfo,$windowSizes)=loadContigBoundFile($contigBoundFile,$windowSize,$windowStep,$verbose) if($contigBoundFile ne "");

my %header2index=();
$header2index{ chrom }=0;
$header2index{ chromStart }=1;
$header2index{ chromEnd }=2;
    
my $chrSort = ($header2index{ chrom }+1).",".($header2index{ chrom }+1);
my $chromStartSort = ($header2index{ chromStart }+1).",".($header2index{ chromStart }+1);
my $sortString="-k".$chrSort." -k".$chromStartSort."n";
    
# sort the bed file

my $inputBedFileName=getFileName($inputBedFile);

print STDERR "checking for sorted file ($inputBedFile) [$sortString]...\n" if($verbose);
my $sortCheck=undef;
if($inputBedFile =~ /\.gz$/) {
    $sortCheck=`gunzip -c '$inputBedFile' | grep -v "track" | sort -c $sortString 2>&1`;
} else {
    $sortCheck=`cat $inputBedFile | grep -v "track" | sort -c $sortString 2>&1`;
}

my $sortedBedFile=$inputBedFile;

chomp($sortCheck);
croak "\nbed file is not sorted!\n\t[$sortCheck]\n" if($sortCheck ne "");

print STDERR "\tsorted!\n" if($verbose);

my $nLines = getNumberOfLines($inputBedFile)-1;
my $progressBucketSize=ceil($nLines / 1000);
my $pcComplete=0;

print STDERR "\n" if($verbose);

my %windowedData=();
my %buffer=();

my $lastChromosome="NA";
my $lastStartPos="NA";
my $lastWindowIndexStart="NA";
my $lastWindowIndexEnd="NA";
my $lineNum=0;

my $regionStart=1;
my $regionEnd=1;

my $negativeFlag=0;

print STDERR "windowing data...\n" if($verbose);
open(IN,inputWrapper($inputBedFile)) or croak "Could not open file [$inputBedFile] - $!";
while(my $line = <IN>) {
    chomp($line);
    next if(($line eq "") or ($line =~ m/^#/));
    
    # skip possible BED headers
    next if($line =~ /^#/);
    next if($line =~ /^track/);
    next if($line =~ /^chrom/);
    
    my @tmp=split(/\t/,$line);
    croak "invalid line # $lineNum - [$line]" if(@tmp < 3);
    
    my $chromosome=$tmp[$header2index{ chrom }];
    my $startPos=$tmp[$header2index{ chromStart }];
    my $endPos=$tmp[$header2index{ chromEnd }];
    
    my $score="NA";
    my $output="NA";
    if(@tmp == 4) {
        $score = $tmp[3] if((defined($tmp[3])) and ($tmp[3] ne "") and ($tmp[3] =~ /^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/));
    }
    if(@tmp > 4) {
        $output = $tmp[3] if((defined($tmp[3])) and ($tmp[3] ne ""));
        $score = $tmp[4] if((defined($tmp[4])) and ($tmp[4] ne "") and ($tmp[4] =~ /^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/));
    }
    
    $regionStart=$startPos if($chromosome ne $lastChromosome);
    $regionEnd=$endPos;
    
    $contigInfo->{$chromosome}->{ start }=$regionStart if( exists($contigInfo->{$chromosome}->{ start }) and ($regionStart < $contigInfo->{$chromosome}->{ start }));
    $contigInfo->{$chromosome}->{ end }=$regionEnd if( exists($contigInfo->{$chromosome}->{ end }) and ($regionEnd > $contigInfo->{$chromosome}->{ end }));
    $contigInfo->{$chromosome}->{ start }=$regionStart if( !exists($contigInfo->{$chromosome}->{ start }) );
    $contigInfo->{$chromosome}->{ end }=$regionEnd if( !exists($contigInfo->{$chromosome}->{ end }) );
    
    $regionStart=$contigInfo->{$chromosome}->{ start } if(exists($contigInfo->{$chromosome}->{ start }));
    $regionEnd=$contigInfo->{$chromosome}->{ end } if(exists($contigInfo->{$chromosome}->{ end }));
    
    croak "(chr) - BED is not properly sorted! [$sortString]" if(($lastChromosome ne "NA") and ($chromosome lt $lastChromosome));
    croak "(start) - BED is not properly sorted! [$sortString]" if(($lastStartPos ne "NA") and ($chromosome eq $lastChromosome) and ($startPos < $lastStartPos));
    
    $startPos -= $regionStart;
    $endPos -= $regionStart;
    
    my $startPos_indexEnd = floor($startPos/($windowSize/$windowStep));
    my $startPos_indexStart = ($startPos_indexEnd-($windowStep-1));
    $startPos_indexStart = 0 if($startPos_indexStart < 0);
    
    my $endPos_indexEnd = floor($endPos/($windowSize/$windowStep));
    my $endPos_indexStart = ($endPos_indexEnd-($windowStep-1));
    $endPos_indexStart = 0 if($endPos_indexStart < 0);
        
    if($lastWindowIndexStart ne "NA") {
        my $bufferBound=$startPos_indexStart;
        $bufferBound=$lastWindowIndexEnd+1 if($chromosome ne $lastChromosome);
        for(my $i=$lastWindowIndexStart;$i<$bufferBound;$i++) {
            my $scoreArrayRef=[];
            $scoreArrayRef=$buffer{$i} if(exists($buffer{$i}));
            
            my $scoreArrayRefStats=listStats($scoreArrayRef) if(@{$scoreArrayRef} > 0);
            my $aggregrateScore="NA";
            $aggregrateScore=$scoreArrayRefStats->{ $windowMode } if(exists($scoreArrayRefStats->{ $windowMode }));
            $windowedData{$lastChromosome}{$i}=$aggregrateScore if($aggregrateScore ne "NA");
            
            $negativeFlag = 1 if(($aggregrateScore ne "NA") and ($aggregrateScore < 0));
            
            delete $buffer{$i} if(exists($buffer{$i}));
        }
    }

    for(my $windowIndex=$startPos_indexStart;$windowIndex<=$endPos_indexEnd;$windowIndex++) {
        push(@{$buffer{$windowIndex}},$score) if($score ne "NA");
        
        $windowSizes->{$chromosome}->{ min } = $windowIndex if(!exists($windowSizes->{$chromosome}->{ min }));
        $windowSizes->{$chromosome}->{ min } = $windowIndex if((exists($windowSizes->{$chromosome}->{ min })) and ($windowIndex < $windowSizes->{$chromosome}->{ min }));
        $windowSizes->{$chromosome}->{ max } = $windowIndex if(!exists($windowSizes->{$chromosome}->{ max }));
        $windowSizes->{$chromosome}->{ max } = $windowIndex if((exists($windowSizes->{$chromosome}->{ max })) and ($windowIndex > $windowSizes->{$chromosome}->{ max }));        
    }
    
    if((($lineNum % $progressBucketSize) == 0) or ($lineNum == ($nLines-1))) {
        $pcComplete = 100 if($lineNum == ($nLines-1));
        print STDERR "\e[A" if(($verbose) and ($lineNum != 0));
        printf STDERR "\t%.2f%% complete ($lineNum/".($nLines-1).")...\n", $pcComplete if($verbose);
        $pcComplete = round((($lineNum/($nLines-1))*100),2);
    }
    
    $lastChromosome=$chromosome;
    $lastStartPos=$startPos;
    $lastWindowIndexStart=$startPos_indexStart;
    $lastWindowIndexEnd=$endPos_indexEnd;
    $lineNum++;    
    
}
close(IN);

# empty last of buffer
if($lastWindowIndexStart ne "NA") {
    for(my $i=$lastWindowIndexStart;$i<($lastWindowIndexEnd+1);$i++) {
        my $scoreArrayRef=[];
        $scoreArrayRef=$buffer{$i} if(exists($buffer{$i}));
        
        my $scoreArrayRefStats=listStats($scoreArrayRef) if(@{$scoreArrayRef} > 0);
        my $aggregrateScore="NA";
        $aggregrateScore=$scoreArrayRefStats->{ $windowMode } if(exists($scoreArrayRefStats->{ $windowMode }));
        $windowedData{$lastChromosome}{$i}=$aggregrateScore if($aggregrateScore ne "NA");
        $negativeFlag = 1 if(($aggregrateScore ne "NA") and ($aggregrateScore < 0));
        
        delete $buffer{$i} if(exists($buffer{$i}));
    }
}
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

print STDERR "writing bed/bedGraph ...\n" if($verbose);

my $fileName=$output."__".$windowSize."__".$windowStep."__".$windowMode;
my $bedFile=$fileName.".bed";
my $bedGraphFile=$fileName.".bedGraph";        
my $trackName = $fileName;

open(BED,outputWrapper($bedFile)) or croak "Could not open file [$bedFile] - $!" if($cisMode == 0);
open(BEDGRAPH,outputWrapper($bedGraphFile)) or croak "Could not open file [$bedGraphFile] - $!" if($cisMode == 0);
print BED "track name='".$trackName."' description='smoothed BED file (".$windowSize."__".$windowStep."__".$windowMode."__".$missingScore.")'\n" if($cisMode == 0);
my $viewOptions="autoScale=on ";
$viewOptions="autoScale=off viewLimits=-".$yBound.":".$yBound if(($cisMode == 0) and (($yBound ne "") and ($negativeFlag)));
$viewOptions="autoScale=off viewLimits=0:".$yBound if(($cisMode == 0) and (($yBound ne "") and ($negativeFlag == 0)));
print BEDGRAPH "track type=bedGraph name='".$trackName."' description='smoothed BED file (".$windowSize."__".$windowStep."__".$windowMode."__".$missingScore.")' visibility=full ".$viewOptions." maxHeightPixels=128:64:32 color=0,0,0 altColor=100,100,100\n" if($cisMode == 0);

foreach my $chromosome ( sort keys %{$windowSizes} ) {
    my $min="NA";
    $min = $windowSizes->{$chromosome}->{ min };
    my $max="NA";
    $max = $windowSizes->{$chromosome}->{ max };
    
    croak "cannot find chromosome bounds - $chromosome" if( (!exists($contigInfo->{$chromosome}->{ start })) or (!exists($contigInfo->{$chromosome}->{ end })) );
    $regionStart=$contigInfo->{$chromosome}->{ start } if(exists($contigInfo->{$chromosome}->{ start }));
    $regionEnd=$contigInfo->{$chromosome}->{ end } if(exists($contigInfo->{$chromosome}->{ end }));
    
    if($cisMode) {
        $fileName=$output."__".$chromosome."__".$windowSize."__".$windowStep."__".$windowMode;
        $trackName = $output."__".$chromosome;
        $bedFile=$fileName.".bed";
        $bedGraphFile=$fileName.".bedGraph";
        
        open(BED,outputWrapper($bedFile)) or croak "Could not open file [$bedFile] - $!";
        open(BEDGRAPH,outputWrapper($bedGraphFile)) or croak "Could not open file [$bedGraphFile] - $!";
        
        print BED "track name='".$trackName."' description='smoothed BED file (".$windowSize."__".$windowStep."__".$windowMode."__".$missingScore.")'\n";
        my $viewOptions="autoScale=on ";
        $viewOptions="autoScale=off viewLimits=-".$yBound.":".$yBound if(($yBound ne "") and ($negativeFlag));
        $viewOptions="autoScale=off viewLimits=0:".$yBound if(($yBound ne "") and ($negativeFlag == 0));
        print BEDGRAPH "track type=bedGraph name='".$trackName."' description='smoothed BED file (".$windowSize."__".$windowStep."__".$windowMode."__".$missingScore.")' visibility=full ".$viewOptions." color=0,0,0 altColor=100,100,100\n";
    }

    for(my $windowIndex=$min;$windowIndex<$max;$windowIndex++) {    
        my $aggregrateScore=$missingScore;
        $aggregrateScore=$windowedData{$chromosome}{$windowIndex} if(exists($windowedData{$chromosome}{$windowIndex}));
        my $windowStart=$regionStart+($windowIndex*$windowStepSize);
        my $windowEnd=($windowStart+$windowSize)-1;
        $windowEnd=$regionEnd if($windowEnd > $regionEnd);
        my $windowName=$chromosome.":".$windowStart."-".$windowEnd;
        
        print BED "$chromosome\t$windowStart\t$windowEnd\t$windowName\t$aggregrateScore\n";
        print BEDGRAPH "$chromosome\t$windowStart\t$windowEnd\t$aggregrateScore\n";
    }
    
    close(BED) if($cisMode);
    close(BEDGRAPH) if($cisMode);

    if($cisMode) {
        my $plotTitle = $fileName;
        $plotTitle =~ s/(...\n{1,30})/$1/gs;

        my $plotBedScriptPath=$scriptPath."plotBED.R";
        system("Rscript '".$plotBedScriptPath."' '".$cwd."' '".$bedFile."' '".$plotTitle."' ".$yBound." > /dev/null");
    }
    
}

print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);