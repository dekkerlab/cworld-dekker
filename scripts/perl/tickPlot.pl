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

    my ($inputMatrix,$verbose,$output,$elementBedFile,$elementZoneSize,$bucketSpan,$excludeCis,$excludeTrans,$excludeZero);
    
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
    
    if( exists($opts->{ elementBedFile }) ) {
        $elementBedFile = $opts->{ elementBedFile };
    } else {
        print STDERR "\nERROR: Option elementBedFile|ebf is required.\n";
        help();
    }
    
    if( exists($opts->{ elementZoneSize }) ) {
        $elementZoneSize = $opts->{ elementZoneSize };
    } else {
        $elementZoneSize=100000;
    }
    
    if( exists($opts->{ bucketSpan }) ) {
        $bucketSpan = $opts->{ bucketSpan };
    } else {
        $bucketSpan=1;
    }
    
    if( exists($opts->{ excludeCis }) ) {
        $excludeCis = 1;
    } else {
        $excludeCis = 0;
    }
    
    if( exists($opts->{ excludeTrans }) ) {
        $excludeTrans = 1;
    } else {
        $excludeTrans = 0;
    }
    
    if( exists($opts->{ excludeZero }) ) {
        $excludeZero = 1;
    } else {
        $excludeZero = 0;
    }
    
    return($inputMatrix,$verbose,$output,$elementBedFile,$elementZoneSize,$bucketSpan,$excludeCis,$excludeTrans,$excludeZero);
}

sub intro() {
    print STDERR "\n";
    
    print STDERR "Tool:\t\telementPileUp.pl\n";
    print STDERR "Version:\t".$cworld::dekker::VERSION."\n";
    print STDERR "Summary:\tpile up cData around specified list of 'elements'\n";
    
    print STDERR "\n";
}

sub help() {
    intro();
    
    print STDERR "Usage: perl elementPileUp.pl [OPTIONS] -i <inputMatrix>\n";
    
    print STDERR "\n";
    
    print STDERR "Required:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-i", "[]", "input matrix file");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ebf", "[]", "element bed file");
    
    print STDERR "\n";
    
    print STDERR "Options:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-v", "[]", "FLAG, verbose mode");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-o", "[]", "prefix for output file(s)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-zs", "[]", "elementZoneSize, size of zone to use around element (x axis - in BP)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--bs", "[1]", "bucketSpan, bucket size for discrete binning");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ec", "[]", "FLAG, exclude CIS data");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--et", "[]", "FLAG, exclude TRANS data");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ez", "[]", "FLAG, ignore zeros from all calculations");    
    
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

sub tickPlot($$$$) {
    my $matrixObject=shift;
    my $pairwiseFile=shift;
    my $elements=shift;
    my $bucketSpan=shift;
    
    my $inc2header=$matrixObject->{ inc2header };
    my $header2inc=$matrixObject->{ header2inc };
    my $numXHeaders=$matrixObject->{ numXHeaders };
    my $numYHeaders=$matrixObject->{ numYHeaders };
    my $numTotalHeaders=$matrixObject->{ numTotalHeaders };
    my $headerSizing=$matrixObject->{ headerSizing };
    my $headerSpacing=$matrixObject->{ headerSpacing };
    my $output=$matrixObject->{ output };
    
    my %elementHeaders=();
    my $numElements=@{$elements};
    for(my $e=0;$e<$numElements;$e++) {
        my $elementHeader=$elements->[$e]->{ name };
        my $elementObject=getHeaderObject($elementHeader,1);
        
        my $elementIndex=-1;
        $elementIndex=$header2inc->{ xy }->{$elementHeader} if(exists($header2inc->{ xy }->{$elementHeader}));
        croak "invalid header ($elementHeader)" if($elementIndex == -1);
        
        $elementHeaders{$elementHeader}=$e;
    }
    
    
    my $lineNum=0;
    my @tmpList=();
    open(IN,inputWrapper($pairwiseFile)) or croak "Could not open file [$pairwiseFile] - $!";
    while(my $line = <IN>) {
        chomp($line);
        next if(($line eq "") or ($line =~ m/^#/));
        
        $lineNum++;
        
        next if($lineNum == 1);
        
        my ($header_1,$header_2,$score)=split(/\t/,$line);
        push(@tmpList,$score) if($score ne "NA");
    }
    close(IN);

    my %tickData=();
    
    #trim off top/bottom % of outliers
    my $tmpListStats=listStats(\@tmpList,0.005);
    undef(@tmpList);
    my $min=$tmpListStats->{ min };
    my $max=$tmpListStats->{ max };
    my $trimmedMin=$tmpListStats->{ trimmedMin };
    my $trimmedMax=$tmpListStats->{ trimmedMax };
    
    $trimmedMin=roundNearest($trimmedMin,$bucketSpan);
    $trimmedMax=roundNearest($trimmedMax,$bucketSpan);
    my $nBuckets=ceil(($trimmedMax-$trimmedMin)/$bucketSpan);
    
    # init tick data
    for(my $ec=-1;$ec<=2;$ec++) {
        for(my $i=0;$i<$nBuckets;$i++) {
            $tickData{$i}{$ec}=0;
        }
    }

    $lineNum=0;
    open(IN,inputWrapper($pairwiseFile)) or croak "Could not open file [$pairwiseFile] - $!";
    while(my $line = <IN>) {
        chomp($line);
        next if(($line eq "") or ($line =~ m/^#/));
        
        $lineNum++;
        
        next if($lineNum == 1);
        
        my ($header_1,$header_2,$score)=split(/\t/,$line);
        
        my $header_1_flag=0;
        $header_1_flag = 1 if(exists($elementHeaders{$header_1}));
        
        my $header_2_flag=0;
        $header_2_flag = 1 if(exists($elementHeaders{$header_2}));
        
        my $bucketIndex=round(($score-$trimmedMin)/$bucketSpan);
        
        $bucketIndex=0 if($bucketIndex < 0);
        $bucketIndex=$nBuckets if($bucketIndex > $nBuckets);
        
        my $elementCount=($header_1_flag+$header_2_flag);
        $tickData{$bucketIndex}{ -1 }++;
        $tickData{$bucketIndex}{$elementCount}++;
        
    }
    close(IN);
    
    my $tickMatrixFile=$output.".tick.matrix.gz";
    open(OUT,outputWrapper($tickMatrixFile)) or croak "Could not open file [$tickMatrixFile] - $!";
    
    for(my $i=0;$i<$nBuckets;$i++) {
        my $bucket=($trimmedMin+($i*$bucketSpan));
        my $bucketStart=$bucket-($bucketSpan/2);
        my $bucketEnd=$bucket+($bucketSpan/2);
        
        $bucketStart="+".$bucketStart if($bucketStart > 0);
        $bucketEnd="+".$bucketEnd if($bucketEnd > 0);
        
        print OUT "\t".$bucketStart."___".$bucketEnd;
    }
    print OUT "\n";
    
    for(my $ec=0;$ec<=2;$ec++) {
        for(my $i=0;$i<$nBuckets;$i++) {
            my $bucket=($trimmedMin+($i*$bucketSpan));
            
            my $allCount=$tickData{$i}{-1};
            my $tmpCount=$tickData{$i}{$ec};
            
            my $tmpPercent=0;
            $tmpPercent = round((($tmpCount/$allCount)*100),2) if($allCount != 0);
            
            print OUT "elementCount_".$ec if($i == 0);
            print OUT "\t".$tmpPercent;
        }
        
        print OUT "\n";
    }
    
    close(OUT);
    
    return($tickMatrixFile);
    
}

    
my %options;
my $results = GetOptions( \%options,'inputMatrix|i=s','verbose|v','output|o=s','elementBedFile|ebf=s','elementZoneSize|ezs=s','bucketSpan|bs=s','excludeCis|ec','excludeTrans|et','excludeZero|ez') or croak help();

my ($inputMatrix,$verbose,$output,$elementBedFile,$elementZoneSize,$bucketSpan,$excludeCis,$excludeTrans,$excludeZero)=check_options( \%options );

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

croak "matrix must be symmetrical!" if($symmetrical == 0);

my $elementFileName=getFileName($elementBedFile);
print STDERR "validating $elementFileName ...\n" if($verbose);
validateBED($elementBedFile);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

# add suffix
$output .= "___".$elementFileName;
$matrixObject->{ output }=$output;

print STDERR "running headers2bed ...\n" if($verbose);
my $headerBEDFile=headers2bed($matrixObject);
print STDERR "\t$headerBEDFile\n" if($verbose);

print STDERR "\n" if($verbose);

print STDERR "intersecting BED files ...\n" if($verbose);
my $bedOverlapFile=intersectBED($headerBEDFile,$elementBedFile,$elementZoneSize);
print STDERR "\t$bedOverlapFile\n" if($verbose);
system("rm '".$headerBEDFile."'");

print STDERR "\n" if($verbose);

print STDERR "loading BED file ...\n" if($verbose);
my ($elements)=loadBED($bedOverlapFile);
print STDERR "\tfound ".@{$elements}." elements\n" if($verbose);
system("rm '".$bedOverlapFile."'");

croak "found no overlapping headers!" if(@{$elements} == 0);

print STDERR "\n" if($verbose);

my $skipNA=1;
print STDERR "printing pairwise file...\n" if($verbose);
my $pairwiseFile=matrix2pairwise($matrixObject,$inputMatrix,$excludeCis,$excludeTrans,$skipNA,$excludeZero);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

tickPlot($matrixObject,$pairwiseFile,$elements,$bucketSpan);