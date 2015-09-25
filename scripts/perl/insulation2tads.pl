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

    my ($insulationFile,$boundaryFile,$verbose,$output,$minBoundaryStrength,$minTadStrength,$yBound);
    
    my $ret={};
    
    if( exists($opts->{ insulationFile }) ) {
        $insulationFile = $opts->{ insulationFile };
    } else {
        print STDERR "\nERROR: Option insulationFile|i is required.\n";
        help();
    }
    
    if( exists($opts->{ boundaryFile }) ) {
        $boundaryFile = $opts->{ boundaryFile };
    } else {
        print STDERR "\nERROR: Option boundaryFile|b is required.\n";
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
    
    if( exists($opts->{ minBoundaryStrength }) ) {
        $minBoundaryStrength = $opts->{ minBoundaryStrength };
    } else {
        $minBoundaryStrength = 0;
    }
    
    if( exists($opts->{ minTadStrength }) ) {
        $minTadStrength = $opts->{ minTadStrength };
    } else {
        $minTadStrength = 0;
    }
    
    if( exists($opts->{ yBound }) ) {
        $yBound = $opts->{ yBound };
    } else {
        $yBound = 0;
    }
    
    $ret->{ insulationFile }=$insulationFile;
    $ret->{ boundaryFile }=$boundaryFile;
    $ret->{ verbose }=$verbose;
    $ret->{ output }=$output;
    $ret->{ minBoundaryStrength }=$minBoundaryStrength;
    $ret->{ minTadStrength }=$minTadStrength;
    $ret->{ yBound }=$yBound;
    
    return($ret,$insulationFile,$boundaryFile,$verbose,$output,$minBoundaryStrength,$minTadStrength,$yBound);
}

sub intro() {
    print STDERR "\n";
    
    print STDERR "Tool:\t\t".$tool."\n";
    print STDERR "Version:\t".$cworld::dekker::VERSION."\n";
    print STDERR "Summary:\tcreate tad specific headers\n";
    
    print STDERR "\n";
}

sub help() {
    intro();
    
    print STDERR "Usage: perl insulation2tads.pl [OPTIONS] -i <inputMatrix> -b <boundaryFile>\n";
    
    print STDERR "\n";
    
    print STDERR "Required:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-i", "[]", "cWorld insulation file");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-b", "[]", "cWorld boundary file");
    
    print STDERR "\n";
    
    print STDERR "Options:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-v", "[]", "FLAG, verbose mode");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-o", "[]", "prefix for output file(s)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--mbs", "[]", "min boundary strength (0-inf)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--mts", "[]", "min tad strength (0-inf)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--yb", "[auto]", "yBound, 0 - yBound for bedGraph");
    
    print STDERR "\n";
    
    print STDERR "Notes:";
    print STDERR "
    This script can assemble consecutive non-overlapping TADs from a insulation/boundary file.\n";
    
    
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

sub loadBoundaries($$;$) {
    # required
    my $boundaryFile=shift;
    my $minBoundaryStrength=shift;
    # optional
    my $verbose=0;
    $verbose=shift if @_;
    
    croak "boundaryFile [$boundaryFile] does not exist" if(!(-e $boundaryFile));
    
    my $groupIndex=0;
    my $lastChromosome="NA";
    
    my $lineNum=0;
    my %header2index=();
    
    my @boundaryStrength_arr=();
    my %boundaries=();
    open(IN,inputWrapper($boundaryFile)) or croak "Could not open file [$boundaryFile] - $!";
    while(my $line = <IN>) {
        chomp($line);
        next if(($line eq "") or ($line =~ m/^#/));
        
        # header  start   end     binStart        binEnd  binMidpoint     boundaryHeader  boundaryStrength
        if($lineNum == 0) {
            my @headers=split(/\t/,$line);
            for(my $i=0;$i<@headers;$i++) {
                my $header=$headers[$i];
                $header="boundaryStrength" if($header eq "insulationScore"); # old version
                $header2index{$header}=$i;
            }
            croak "invalid boundary file! (incorrect column headers!) [missing 'header']" if(!exists($header2index{ header }));
            croak "invalid boundary file! (incorrect column headers!) [missing 'start']" if(!exists($header2index{ start }));
            croak "invalid boundary file! (incorrect column headers!) [missing 'end']" if(!exists($header2index{ end }));
            croak "invalid boundary file! (incorrect column headers!) [missing 'boundaryStrength']" if(!exists($header2index{ boundaryStrength }));
            $lineNum++;
            next;
        }
        
        my @tmp=split(/\t/,$line);
        my $boundaryHeader=$tmp[ $header2index{ header } ];
        my $boundaryStart=$tmp[ $header2index{ start } ];
        my $boundaryEnd=$tmp[ $header2index{ end } ];
        my $boundaryStrength=$tmp[ $header2index{ boundaryStrength } ];
        
        next if(($boundaryStrength ne "NA") and ($boundaryStrength < $minBoundaryStrength));
        push(@boundaryStrength_arr,$boundaryStrength) if($boundaryStrength ne "NA");
        
        my $boundaryHeaderInfo=getHeaderObject($boundaryHeader);
        my $boundaryHeaderChromosome=$boundaryHeaderInfo->{ chromosome };
        
        $groupIndex = 0 if($lastChromosome ne $boundaryHeaderChromosome);
        
        $boundaries{$boundaryHeaderChromosome}{$groupIndex}{ start }=$boundaryStart;
        $boundaries{$boundaryHeaderChromosome}{$groupIndex}{ end }=$boundaryEnd;
        $boundaries{$boundaryHeaderChromosome}{$groupIndex}{ strength }=$boundaryStrength;
        
        $groupIndex++;
        $lastChromosome=$boundaryHeaderChromosome;
        
    }
    close(IN);
    
    my $boundaryStrengthsStats=listStats(\@boundaryStrength_arr) if(@boundaryStrength_arr > 0);
    my $boundaryStrength_mean="NA";
    $boundaryStrength_mean=$boundaryStrengthsStats->{ median } if(exists($boundaryStrengthsStats->{ median }));
    
    print STDERR "\tnum boundaries = ".scalar @boundaryStrength_arr."\n" if($verbose);
    print STDERR "\tmedian strength = $boundaryStrength_mean\n" if($verbose);
    
    
    return(\%boundaries);
}

sub assembleTADs($$$;$) {
    # required
    my $insulationFile=shift;
    my $boundaries=shift;
    my $inputName=shift;
    # optional
    my $verbose=0;
    $verbose=shift if @_;
    
    croak "insulationFile [$insulationFile] does not exist" if(!(-e $insulationFile));
    
    my %header2group=();
    my %header2groupIndex=();
    my %header2groupOffset=();
    my %groupInfo=();
    my %usableTADs=();
    
    my $NAcount=0;
    my $tadSize=0;
    my $tadStart=-1;
    my $tadEnd=-1;
    
    my $groupOffset=0;
    my $lastGroupOffset=0;
    my $groupIndex=0;
    my $lastBoundaryIndex=-1;
    my $lastChromosome="NA";
    my $assembly="NA";
    
    my $lineNum=0;
    my %header2index=();
    my @tmpTadStrengths=();
    open(IN,inputWrapper($insulationFile)) or croak "Could not open file [$insulationFile] - $!";
    while(my $line = <IN>) {
        chomp($line);
        next if(($line eq "") or ($line =~ m/^#/));
        
        #header  start   end     midpoint        binStart        binEnd  binMidpoint     rawInsulationScore      insulationScore delta   deltaSquare
        if($lineNum == 0) {
            my @headers=split(/\t/,$line);
            for(my $i=0;$i<@headers;$i++) {
                my $header=$headers[$i];
                $header2index{$header}=$i;
            }
            croak "invalid insulation file! (missing header/insulationScore columns" if( (!exists($header2index{ header })) or (!exists($header2index{ insulationScore })) );

            $lineNum++;
            next;
        }
        
        my @tmp=split(/\t/,$line);
        my $header = $tmp[ $header2index{ header } ];
        my $insulationScore = $tmp [ $header2index{ insulationScore } ];
        
        my $headerObject=getHeaderObject($header);
        my $headerSubName=$headerObject->{ subName };
        $assembly=$headerObject->{ assembly };        
        my $headerChromosome=$headerObject->{ chromosome };
        my $headerCoordinates=$headerObject->{ coordinates };
        my $headerStart=$headerObject->{ start };
        my $headerEnd=$headerObject->{ end };
        my $headerMid=round(($headerStart+$headerEnd)/2);
        
        # init tad coordinates
        $tadStart=$headerStart if($tadStart == -1);
        $tadEnd=$headerEnd if($tadEnd == -1);
        
        #print STDERR "$header\t$groupIndex\t$tadStart\t$tadEnd\t$insulationScore\t$NAcount\t$tadSize\n";
        
        # get the boundaries for chromosome
        my $boundarySubHash=$boundaries->{$headerChromosome};
        my $nBoundaries=keys %{$boundarySubHash};
        
        # reset boundary offset when new chr
        $groupOffset = 0 if($lastChromosome ne $headerChromosome);
        $lastGroupOffset = 0 if($lastChromosome ne $headerChromosome);
        
        my $boundaryFlag=0;
        # search for a boundary overlap
        for(my $i=$groupOffset;$i<$nBoundaries;$i++) {
            my $boundaryStart=$boundarySubHash->{$i}->{ start };
            my $boundaryEnd=$boundarySubHash->{$i}->{ end };
            my $boundaryStrength=$boundarySubHash->{$i}->{ strength };
                        
            $groupOffset++ if($boundaryEnd <= $headerStart); # skip boundaries on left
            next if($boundaryEnd <= $headerStart); # skip boundaries on left
            last if($boundaryStart >= $headerEnd); # skip boundaries on right
            
            # found a boundary!
            
            my $tmpTadStrengthStats=listStats(\@tmpTadStrengths) if(@tmpTadStrengths > 0);
            my $tmpTadStrength_min="NA";
            $tmpTadStrength_min=$tmpTadStrengthStats->{ min } if(exists($tmpTadStrengthStats->{ min }));
            my $tmpTadStrength_max="NA";
            $tmpTadStrength_max=$tmpTadStrengthStats->{ max } if(exists($tmpTadStrengthStats->{ max }));
            
            my $leftInsulation="NA";
            $leftInsulation=$tmpTadStrengths[0] if(@tmpTadStrengths > 0);
            my $rightInsulation="NA";
            $rightInsulation=$tmpTadStrengths[-1] if(@tmpTadStrengths > 0);
            my $tmpTadStrength="NA";
            $tmpTadStrength=min(($tmpTadStrength_max-$leftInsulation),($tmpTadStrength_max-$rightInsulation)) if(($leftInsulation ne "NA") or ($rightInsulation ne "NA"));
            
            $groupInfo{$groupIndex}{ size } = $tadSize;
            $groupInfo{$groupIndex}{ NAcount } = $NAcount;
            $groupInfo{$groupIndex}{ strength } = $tmpTadStrength;
            $groupInfo{$groupIndex}{ assembly } = $assembly;
            $groupInfo{$groupIndex}{ chromosome } = $lastChromosome;
            $groupInfo{$groupIndex}{ start } = $tadStart;
            $groupInfo{$groupIndex}{ end } = $tadEnd;
            
            my $usableTADFlag=1;
            $usableTADFlag=0 if($NAcount > ($tadSize/4));
            $usableTADs{$groupIndex}=$usableTADFlag;
            
            $groupOffset=$i;
            $boundaryFlag=1;
            $groupIndex++ if($groupOffset != $lastGroupOffset);
            $lastGroupOffset=$i;
        }
        
        $NAcount++ if($insulationScore eq "NA");
        $tadSize++;
        push(@tmpTadStrengths,$insulationScore) if($insulationScore ne "NA");
        $tadEnd=$headerEnd;
        
        # append on header prefix
        my $groupName="tad".$groupIndex;
        $groupName="boundary".$groupIndex if($boundaryFlag);
        
        $NAcount=0 if($boundaryFlag);
        $tadSize=0 if($boundaryFlag);
        @tmpTadStrengths=() if($boundaryFlag);
        $tadStart=-1 if($boundaryFlag);
        $tadEnd=-1 if($boundaryFlag);
        
        $header2group{$header}=$groupName;
        $header2groupIndex{$header}=$groupIndex;
        $header2groupOffset{$header}=$groupOffset;
        $lastChromosome=$headerChromosome;
        
        $lineNum++;

    }
    close(IN);
    
    my $tmpTadStrengthStats=listStats(\@tmpTadStrengths) if(@tmpTadStrengths > 0);
    my $tmpTadStrength_min="NA";
    $tmpTadStrength_min=$tmpTadStrengthStats->{ min } if(exists($tmpTadStrengthStats->{ min }));
    my $tmpTadStrength_max="NA";
    $tmpTadStrength_max=$tmpTadStrengthStats->{ max } if(exists($tmpTadStrengthStats->{ max }));
    my $tmpTadStrength="NA";
    $tmpTadStrength=($tmpTadStrength_max-$tmpTadStrength_min) if(($tmpTadStrength_max ne "NA") and ($tmpTadStrength_min ne "NA"));

    $groupInfo{$groupIndex}{ size } = $tadSize;
    $groupInfo{$groupIndex}{ NAcount } = $NAcount;
    $groupInfo{$groupIndex}{ strength } = $tmpTadStrength;
    $groupInfo{$groupIndex}{ assembly } = $assembly;
    $groupInfo{$groupIndex}{ chromosome } = $lastChromosome;
    $groupInfo{$groupIndex}{ start } = $tadStart;
    $groupInfo{$groupIndex}{ end } = $tadEnd;
        
    my $usableTADFlag=1;
    $usableTADFlag=0 if($NAcount > ($tadSize/4));
    $usableTADs{$groupIndex}=$usableTADFlag;
    
    
    
    my $headerMapFile=$inputName.".tadHeaders.map";
    open(MAP,outputWrapper($headerMapFile)) or croak "Could not open file [$headerMapFile] - $!";
    
    $lastChromosome="NA";
    $lineNum=0;
    open(IN,inputWrapper($insulationFile)) or croak "Could not open file [$insulationFile] - $!";
    while(my $line = <IN>) {
        chomp($line);
        next if(($line eq "") or ($line =~ m/^#/));
        
        #header  start   end     midpoint        binStart        binEnd  binMidpoint     rawInsulationScore      insulationScore delta   deltaSquare
        if($lineNum == 0) {
            my @headers=split(/\t/,$line);
            for(my $i=0;$i<@headers;$i++) {
                my $header=$headers[$i];
                $header2index{$header}=$i;
            }
            croak "invalid insulation file! (missing header/insulationScore columns)" if( (!exists($header2index{ header })) or (!exists($header2index{ insulationScore })) );

            $lineNum++;
            next;
        }
        
        my @tmp=split(/\t/,$line);
        my $header = $tmp[ $header2index{ header } ];
        my $insulationScore = $tmp [ $header2index{ insulationScore } ];
        
        my $headerObject=getHeaderObject($header);
        my $headerSubName=$headerObject->{ subName };
        $assembly=$headerObject->{ assembly };        
        my $headerChromosome=$headerObject->{ chromosome };
        my $headerCoordinates=$headerObject->{ coordinates };
        my $headerStart=$headerObject->{ start };
        my $headerEnd=$headerObject->{ end };
        my $headerMid=round(($headerStart+$headerEnd)/2);
        
        my $groupName="";
        $groupName="-".$header2group{$header} if(exists($header2group{$header}));
        my $groupIndex="NA";
        $groupIndex=$header2groupIndex{$header} if(exists($header2groupIndex{$header}));
        my $groupOffset="NA";
        $groupOffset=$header2groupOffset{$header} if(exists($header2groupOffset{$header}));
        
        croak "header ($header) does not exist in boundary file!" if(($groupIndex eq "NA") or ($groupOffset eq "NA"));
        
        $groupName="" if(($usableTADs{$groupIndex} == 0) and ($groupName =~ /^-tad/));
        my $usableTADFlag=0;
        $usableTADFlag=$usableTADs{$groupIndex} if(exists($usableTADs{$groupIndex}));
        
        my $newHeader=$headerSubName."|".$assembly."|".$headerChromosome.$groupName.":".$headerStart."-".$headerEnd;
        print MAP "$header\t$newHeader\n";
        
        $lastBoundaryIndex=$groupIndex;
        $lastChromosome=$headerChromosome;
        
    }
    close(IN);
    
    close(MAP);
    
    return($headerMapFile,\%groupInfo,\%usableTADs);
    
}

sub outputTADs($$$$$;$) {
    # required
    my $inputName=shift;
    my $groupInfo=shift;
    my $usableTADs=shift;
    my $yBound=shift;
    my $minTadStrength=shift;
    # optional
    my $verbose=0;
    $verbose=shift if @_;
    
    # now output tad information
    
    my $tadBedFile=$inputName.".tads.bed";
    open(BED,outputWrapper($tadBedFile)) or croak "Could not open file [$tadBedFile] - $!";
    print BED "track name='".$inputName."-bed' description='".$inputName." bed tads' visibility=squish\n";
    
    my $tadBedGraphFile=$inputName.".tads.bedGraph";
    open(BEDGRAPH,outputWrapper($tadBedGraphFile)) or croak "Could not open file [$tadBedGraphFile] - $!";
    print BEDGRAPH "track type=bedGraph name='".$inputName."-bedGraph' description='".$inputName." bedGraph tads' visibility=dense autoScale=off viewLimits=0:".$yBound." color=0,0,0 altColor=100,100,100\n" if($yBound != 0);
    print BEDGRAPH "track type=bedGraph name='".$inputName."-bedGraph' description='".$inputName." bedGraph tads' visibility=dense autoScale=on color=0,0,0 altColor=100,100,100\n" if($yBound == 0);
    
    my $nTads=0;
    my $nFilteredTads=0;
    foreach my $groupIndex (sort { $a <=> $b } keys(%{$groupInfo}) ) {
        
        my $tadSize=$groupInfo->{$groupIndex}->{ size };
        my $tadStrength=$groupInfo->{$groupIndex}->{ strength };
        my $tadAssembly=$groupInfo->{$groupIndex}->{ assembly };
        my $tadChromosome=$groupInfo->{$groupIndex}->{ chromosome };
        my $tadStart=$groupInfo->{$groupIndex}->{ start };
        my $tadEnd=$groupInfo->{$groupIndex}->{ end };
        my $tadName="tad".$groupIndex."|".$tadAssembly."|".$tadChromosome.":".$tadStart."-".$tadEnd;
        $tadStrength="NaN" if($tadStrength eq "NA");
        
        # strip off chr group if exists for proper UCSC usage
        $tadChromosome=stripChromosomeGroup($tadChromosome);
        
        $nTads++;
        next if($tadStrength < $minTadStrength);
        $nFilteredTads++;
        
        print BED "$tadChromosome\t$tadStart\t$tadEnd\t$tadName\t1000\n" if($usableTADs->{$groupIndex});
        print BEDGRAPH "$tadChromosome\t$tadStart\t$tadEnd\t$tadStrength\n" if($usableTADs->{$groupIndex});
    }
    
    print STDERR "\t$nTads tads\n" if($verbose);
    print STDERR "\t$nFilteredTads filtered tads\n" if($verbose);
    
    close(MAP);
    close(BED);
    
    return($tadBedFile,$tadBedGraphFile);
}
    

my %options;
my $results = GetOptions( \%options,'insulationFile|i=s','boundaryFile|b=s','verbose|v','output|o=s','minBoundaryStrength|mbs=s','minTadStrength|mts=s','yBound|yb=f') or croak help();
my ($ret,$insulationFile,$boundaryFile,$verbose,$output,$minBoundaryStrength,$minTadStrength,$yBound)=check_options( \%options );

intro() if($verbose);

#get the absolute path of the script being used
my $cwd = getcwd();
my $fullScriptPath=abs_path($0);
my @fullScriptPathArr=split(/\//,$fullScriptPath);
@fullScriptPathArr=@fullScriptPathArr[0..@fullScriptPathArr-3];
my $scriptPath=join("/",@fullScriptPathArr);

croak "insulationFile [$insulationFile] does not exist" if(!(-e $insulationFile));
croak "boundaryFile [$boundaryFile] does not exist" if(!(-e $boundaryFile));

$output=getFileName($insulationFile) if($output eq "");

# add suffix
$output .= "--mbs".$minBoundaryStrength."--mts".$minTadStrength;

print STDERR "loading boundaries ...\n" if($verbose);
my $boundaries=loadBoundaries($boundaryFile,$minBoundaryStrength,$verbose);

print STDERR "\n" if($verbose);

print STDERR "assembling tads ...\n" if($verbose);
my ($headerMapFile,$groupInfo,$usableTADs)=assembleTADs($insulationFile,$boundaries,$output,$verbose);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

print STDERR "outputing tads...\n" if($verbose);
my ($tadBedFile,$tadBedGraphFile)=outputTADs($output,$groupInfo,$usableTADs,$yBound,$minTadStrength,$verbose);

print STDERR "\n" if($verbose);