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

    my ($inputInsulation_1,$inputInsulation_2,$verbose,$output,$yBound,$transparentBGFlag);
    
    if( exists($opts->{ inputInsulation_1 }) ) {
        $inputInsulation_1 = $opts->{ inputInsulation_1 };
    } else {
        print STDERR "\nERROR: Option inputInsulation_1|1 is required.\n";
        help();
    }
    
    if( exists($opts->{ inputInsulation_2 }) ) {
        $inputInsulation_2 = $opts->{ inputInsulation_2 };
    } else {
        print STDERR "\nERROR: Option inputInsulation_1|1 is required.\n";
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
        $output = "insulation-difference";
    }
    
    if( exists($opts->{ yBound }) ) {
        $yBound = $opts->{ yBound };
    } else {
        $yBound = 0;
    }
    
    if( exists($opts->{ transparentBGFlag }) ) {
        $transparentBGFlag = 1;
    } else {
        $transparentBGFlag = 0;
    }
    
    return($inputInsulation_1,$inputInsulation_2,$verbose,$output,$yBound,$transparentBGFlag);
}

sub intro() {
    print STDERR "\n";
    
    print STDERR "Tool:\t\tcompareInsulation.pl\n";
    print STDERR "Version:\t".$cworld::dekker::VERSION."\n";
    print STDERR "Summary:\tcompare insulation vector - calculate difference\n";
    
    print STDERR "\n";
}

sub help() {
    intro();
    
    print STDERR "Usage: perl compareInsulation.pl [OPTIONS] -i1 <insulation_1> -i2 <insulation_2>\n";
    
    print STDERR "\n";
    
    print STDERR "Required:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-1", "[]", "insulation vectror 1 file");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-2", "[]", "insulation vectror 2 file");
    
    print STDERR "\n";
    
    print STDERR "Options:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-v", "[]", "FLAG, verbose mode");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-o", "[]", "prefix for output file(s)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--yb", "[]", "optional, fix y-limit for insulation image");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--bg", "[]", "optional, enable transparent background for insulation image");

    print STDERR "\n";
    
    print STDERR "Notes:";
    print STDERR "
    This script can compare two insulation vectors [difference].
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

sub outputInsulationDifferential($$$$$$) {
    my $insulationDifferentialRef=shift;
    my $insulationHeaders=shift;
    my $numHeaders=shift;
    my $output=shift;
    my $headerSpacing=shift;
    my $yBound=shift;
    
    my $outputName=getFileName($output);
    
    open(OUT,outputWrapper($output)) or croak "Could not open file [$output] - $!";
    print OUT "header\tstart\tend\tmidpoint\tbinStart\tbinEnd\tbinMidpoint\tinsulationDifferential\n";
    open(BEDGRAPH,outputWrapper($outputName.".bedGraph")) or croak "Could not open file [$output] - $!";
  
    print BEDGRAPH "track type=bedGraph name='".$output."' description='".$output." - insutation score' visibility=full autoScale=off viewLimits=-".$yBound.":".$yBound." color=255,0,0 altColor=0,0,255\n" if($yBound != 0);
    print BEDGRAPH "track type=bedGraph name='".$output."' description='".$output." - insutation score' visibility=full autoScale=on color=255,0,0 altColor=0,0,255\n" if($yBound == 0);
    
    my $lastChromosome="NA";
    for(my $y=0;$y<$numHeaders;$y++) {
        # dump insulation data to file
        my $yHead=$insulationHeaders->{ y }->{$y};
        my $insulationDifferential="NA";
        $insulationDifferential=$insulationDifferentialRef->{$yHead} if(exists($insulationDifferentialRef->{$yHead}));
        
        my $yHeadObject=getHeaderObject($yHead);
        my $yHeadChromosome=$yHeadObject->{ chromosome };
        my $yHeadStart=$yHeadObject->{ start };
        my $yHeadEnd=$yHeadObject->{ end };
        my $yHeadMidpoint=round(($yHeadStart+$yHeadEnd)/2);
        
        my $binStart = round($yHeadStart/$headerSpacing);
        my $binEnd = round($yHeadEnd/$headerSpacing);
        my $binMidpoint=(($binStart+$binEnd)/2);
        
        print OUT "$yHead\t$yHeadStart\t$yHeadEnd\t$yHeadMidpoint\t$binStart\t$binEnd\t$binMidpoint\t$insulationDifferential\n";
        
        $insulationDifferential=0 if($insulationDifferential eq "NA");
        print BEDGRAPH "$yHeadChromosome\t$yHeadStart\t$yHeadEnd\t$insulationDifferential\n";
        $lastChromosome=$yHeadChromosome;
    }
    
    close(OUT);
    close(BEDGRAPH);
}

sub getInsulationData($) {
    # extract headers from an insulation file
    
    my $inputInsulation=shift;
    
    my %insulationHeaders=();
    my %insulationHeadersIndex=();
    my %insulationData=();

    my %header2index=();
    
    my $lineNum=0;
    my $headerInc=0;
    open(IN,inputWrapper($inputInsulation)) or croak "Could not open file [$inputInsulation] - $!";
    while(my $line = <IN>) {
        chomp($line);
        next if(($line eq "") or ($line =~ m/^#/));
        
        my @tmp=split(/\t/,$line);
        
        if($lineNum == 0) {
            for(my $i=0;$i<@tmp;$i++) {
                $header2index{$tmp[$i]}=$i;
                #print STDERR "$i\t$tmp[$i]\n";
            }
        } else {
            my $header=$tmp[$header2index{ header }];
            my $midpoint=$tmp[$header2index{ midpoint }];
            my $score=$tmp[$header2index{ insulationScore }];
            
            my $headerObject=getHeaderObject($header,1);
            my $headerSubName=$headerObject->{ subName };
            my $deGroupedHeader=deGroupHeader($header,"liteChr",$lineNum);
            my $deGroupedHeaderObject=getHeaderObject($deGroupedHeader,1);
            my $Header=$deGroupedHeaderObject->{ coords };
            
            # de group the header (if group exists)
            $header=$deGroupedHeader;

            $insulationHeaders{$header}=$headerInc;
            $insulationHeadersIndex{$headerInc}=$header;
            $insulationData{$headerInc}{ header }=$header;
            $insulationData{$headerInc}{ midpoint }=$midpoint;
            $insulationData{$headerInc}{ score }=$score;
            $headerInc++;
        }
        
        $lineNum++;
    }
    close(IN);
    
    return(\%insulationHeaders,\%insulationHeadersIndex,\%insulationData);
    
}

sub calculateInsulationDifferential($$$) {
    my $insulationData_1=shift;
    my $insulationData_2=shift;
    my $numHeaders=shift;
    
    my %insulationDifferential=();
    
    for(my $i=0;$i<$numHeaders;$i++) {
        
        my $insulationHeader1=$insulationData_1->{$i}->{ header };
        my $insulationHeader2=$insulationData_2->{$i}->{ header };
        
        my $insulationScore1=$insulationData_1->{$i}->{ score };
        my $insulationScore2=$insulationData_2->{$i}->{ score };

        my $insulationDifferential="NA";
        if(($insulationScore1 ne "NA") and ($insulationScore2 ne "NA")) {
            $insulationDifferential=($insulationScore1-$insulationScore2);
        }
        
        $insulationDifferential{$insulationHeader1}=$insulationDifferential;
        
    }

    return(\%insulationDifferential);
}

my %options;
my $results = GetOptions( \%options,'inputInsulation_1|1=s','inputInsulation_2|2=s','verbose|v','output|o=s','yBound|yb=f','transparentBGFlag|bg') or croak help();

my ($inputInsulation_1,$inputInsulation_2,$verbose,$output,$yBound,$transparentBGFlag)=check_options(\%options );

intro() if($verbose);

#get the absolute path of the script being used
my $cwd = getcwd();
my $fullScriptPath=abs_path($0);
my @fullScriptPathArr=split(/\//,$fullScriptPath);
@fullScriptPathArr=@fullScriptPathArr[0..@fullScriptPathArr-3];
my $scriptPath=join("/",@fullScriptPathArr);

croak "inputInsulation_1 [$inputInsulation_1] does not exist" if(!(-e $inputInsulation_1));
croak "inputInsulation_2 [$inputInsulation_2] does not exist" if(!(-e $inputInsulation_2));
croak "inputInsulation_1 == inputInsulation_2" if($inputInsulation_1 eq $inputInsulation_2);

print STDERR "parsing insulation1 headers..\n" if($verbose);
my $insulationHeaders_1={};
my $insulationHeadersIndex_1={};
my $insulationData_1={};
($insulationHeaders_1,$insulationHeadersIndex_1,$insulationData_1)=getInsulationData($inputInsulation_1);
my $numHeaders_1=keys(%{$insulationHeaders_1});
print STDERR "\tfound $numHeaders_1\n" if($verbose);
print STDERR "\n" if($verbose);

print STDERR "parsing insulation2 headers..\n" if($verbose);
my $insulationHeaders_2={};
my $insulationHeadersIndex_2={};
my $insulationData_2={};
($insulationHeaders_2,$insulationHeadersIndex_2,$insulationData_2)=getInsulationData($inputInsulation_2);
my $numHeaders_2=keys(%{$insulationHeaders_2});
print STDERR "\tfound $numHeaders_2\n" if($verbose);

print STDERR "\n" if($verbose);

print STDERR "validating input files...\n" if($verbose);
# lazy test for identical insulation input structure
croak "matrix is not symmetrical" if($numHeaders_1 != $numHeaders_2);
my $numHeaders=$numHeaders_1=$numHeaders_2;

# enforce idential insulation input structure 
foreach my $header ( keys(%{$insulationHeaders_1}) ) {
    if(!exists($insulationHeaders_2->{$header})) {
        croak "input1 ($header) does not exist in input2-headers."
    }
}
foreach my $header ( keys(%{$insulationHeaders_2}) ) {
    if(!exists($insulationHeaders_1->{$header})) {
        croak "input2 ($header) does not exist in input1-headers."
    }
}

# headers are the same, use one going foward
my $insulationHeaders={};
$insulationHeaders->{ y }=$insulationHeaders_1=$insulationHeaders_2;
$insulationHeaders->{ x }=$insulationHeaders_1=$insulationHeaders_2;
my $insulationHeadersIndex={};
$insulationHeadersIndex->{ y }=$insulationHeadersIndex_1=$insulationHeadersIndex_2;
$insulationHeadersIndex->{ x }=$insulationHeadersIndex_1=$insulationHeadersIndex_2;

my $imageWidth=$numHeaders*1.1;
$imageWidth=900 if($imageWidth < 900);

print STDERR "\tdone\n" if($verbose);

# assume headers are symmetrical at this point

# get fragment spacing (i.e. bin size)
my ($equalSpacingFlag,$equalSizingFlag,$headerSpacing,$headerSizing)=getHeaderSpacing($insulationHeadersIndex->{ y });

print STDERR "\n" if($verbose);

my $inputInsulationName_1=getFileName($inputInsulation_1);
my $inputInsulationName_2=getFileName($inputInsulation_2);
$output = $inputInsulation_1."___".$inputInsulation_2 if($output eq "");

# add suffix
$output .= ".differential.gz";

# calculate difference between insulation data
print STDERR "calculating insulation differences...\n" if($verbose);
my ($insulationDifferential)=calculateInsulationDifferential($insulationData_1,$insulationData_2,$numHeaders);
print STDERR "\tdone\n" if($verbose);
print STDERR "\n" if($verbose);

# write insulation differential data to file
print STDERR "outputing insulation...\n" if($verbose);
outputInsulationDifferential($insulationDifferential,$insulationHeadersIndex,$numHeaders,$output,$headerSpacing,$yBound);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

# plot insulation using R
print STDERR "plotting insulation...\n" if($verbose);
system("Rscript '".$scriptPath."/R/compareInsulation.R' '".$cwd."' '".$output."' '".$inputInsulationName_1."' '".$inputInsulationName_2."' ".$headerSizing." ".$headerSpacing." ".$imageWidth." ".$yBound." ".$transparentBGFlag." > /dev/null");
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);
