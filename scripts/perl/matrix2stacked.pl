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

    my ($inputMatrix,$verbose,$output,$minDistance,$maxDistance);
    
    my $ret={};
    
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
    
    if( exists($opts->{ minDistance }) ) {
        $minDistance = $opts->{ minDistance };
    } else {
        $minDistance=undef;
    }
    
    if( exists($opts->{ maxDistance }) ) {
        $maxDistance = $opts->{ maxDistance };
    } else {
        $maxDistance = 40000;
    }
    
    $ret->{ inputMatrix }=$inputMatrix;
    $ret->{ verbose }=$verbose;
    $ret->{ output }=$output;
    $ret->{ minDistance }=$minDistance;
    $ret->{ maxDistance }=$maxDistance;
    
    return($ret,$inputMatrix,$verbose,$output,$minDistance,$maxDistance);
}


sub intro() {
    print STDERR "\n";
    
    print STDERR "Tool:\t\t".$tool."\n";
    print STDERR "Version:\t".$cworld::dekker::VERSION."\n";
    print STDERR "Summary:\ttransform matrix into stacked anchor matrix\n";
    
    print STDERR "\n";
}

sub help() {
    intro();
    
    print STDERR "Usage: perl matrix2stacked.pl [OPTIONS] -i <inputMatrix>\n";
    
    print STDERR "\n";
    
    print STDERR "Required:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-i", "[]", "input matrix file (can accept multiple files)");
    
    print STDERR "\n";
    
    print STDERR "Options:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-v", "[]", "FLAG, verbose mode");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-o", "[]", "prefix for output file(s)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--minDist", "[]", "minimum allowed interaction distance, exclude < N distance (in BP)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--maxDist", "[]", "maximum allowed interaction distance, exclude > N distance (in BP)");
    
    print STDERR "\n";
    
    print STDERR "Notes:";
    print STDERR "
    This script transforms a matrix into stacked anchor matrix.
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

sub subsetData($$$$) {
    # required
    my $inputMatrix=shift;
    my $header2inc=shift;
    my $minDistance=shift;
    my $maxDistance=shift;
    
    my $matrixObject=getMatrixObject($inputMatrix);
    my $symmetrical=$matrixObject->{ symmetrical };
    
    my $lineNum=0;
    my @xHeaders=();
    
    my %matrix=();
    my $maxLeftSubsetSize=0;
    my $maxRightSubsetSize=0;
    
    open(IN,inputWrapper($inputMatrix)) or croak "Could not open file [$inputMatrix] - $!";
    while(my $line = <IN>) {
        chomp($line);
        next if(($line eq "") or ($line =~ m/^#/));
        
        my $leftSubsetSize=0;
        my $rightSubsetSize=0;
        
        if($lineNum == 0) {
            @xHeaders=split(/\t/,$line);
        } else {
            my @data=split(/\t/,$line);
            my $dsize=@data;
            my $yHeader=$data[0];
            my $yHeaderObject=getHeaderObject($yHeader);
            my $yHeaderMidpoint=$yHeaderObject->{ midpoint };
            
            my $yIndex=-1;
            $yIndex = $header2inc->{ y }->{$yHeader} if(defined($header2inc->{ y }->{$yHeader}));
            print STDERR "\nWARNING - header ($yHeader) does not exists in header2inc!\n\n" if($yIndex == -1);
            next if($yIndex == -1);
                
            for(my $d=1;$d<$dsize;$d++) {
                my $xHeader=$xHeaders[$d];
                my $xHeaderObject=getHeaderObject($xHeader);
                my $xHeaderMidpoint=$xHeaderObject->{ midpoint };
                
                my $cScore=$data[$d];
                # skip if cScore is not a valid number
                $cScore = "NA" if($cScore !~ (/^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/));
                $cScore = sprintf("%.4f",$cScore) if($cScore ne "NA");
                
                my $xIndex=-1;
                $xIndex = $header2inc->{ x }->{$xHeader} if(defined($header2inc->{ x }->{$xHeader}));
                print STDERR "\nWARNING - header ($xHeader) does not exists in header2inc!\n\n" if($xIndex == -1);
                next if($xIndex == -1);
                
                my $interactionDistance=getInteractionDistance($matrixObject,$yHeaderObject,$xHeaderObject,1);
                next if($interactionDistance == -1);
                next if($interactionDistance > $maxDistance);
                $interactionDistance *= -1 if($xHeaderMidpoint < $yHeaderMidpoint);
                
                $leftSubsetSize++ if($interactionDistance < 0);
                $rightSubsetSize++ if($interactionDistance > 0);
                
                # ensure symmetrical data
                if(($symmetrical) and (exists($matrix{$xIndex}{$yIndex}))) {
                    print STDERR "\nERROR - data is not symmetrical ($xIndex,$yIndex) [$cScore vs ".$matrix{$xIndex}{$yIndex} ."]\n\n" if($matrix{$xIndex}{$yIndex} ne $cScore);
                }
                    
                $matrix{$yIndex}{$xIndex}=$cScore;
            }
            $maxLeftSubsetSize=$leftSubsetSize if($leftSubsetSize > $maxLeftSubsetSize);
            $maxRightSubsetSize=$rightSubsetSize if($rightSubsetSize > $maxRightSubsetSize);
        }
        $lineNum++;
    }
    close(IN);

    my $maxSubsetSize=max($maxLeftSubsetSize,$maxRightSubsetSize);
    
    return(\%matrix,$maxSubsetSize);
}

sub matrix2stacked($$$$$$) {
    my $matrixObject=shift;
    my $matrix=shift;
    my $minDistance=shift;
    my $maxDistance=shift;
    my $missingValue=shift;
    my $subsetSpan=shift;
    
    my $inc2header=$matrixObject->{ inc2header };
    my $numYHeaders=$matrixObject->{ numYHeaders };
    my $numXHeaders=$matrixObject->{ numXHeaders };
    my $headerSizing=$matrixObject->{ headerSizing };
    my $headerSpacing=$matrixObject->{ headerSpacing };
    my $verbose=$matrixObject->{ verbose };
    my $output=$matrixObject->{ output };
    
    my $stackedMatrixFile=$output.".stacked.matrix.gz";
    open(OUT,outputWrapper($stackedMatrixFile)) or croak "Could not open file [$stackedMatrixFile] - $!";
    
    my $nOffsetBins=floor(($maxDistance-($headerSizing-$headerSpacing))/$headerSpacing);
    $nOffsetBins=$subsetSpan;
    
    for(my $i=-$nOffsetBins;$i<=$nOffsetBins;$i++) {
        my $binOffset=(($i+($headerSizing-$headerSpacing))*$headerSpacing);
        print OUT "\t$binOffset";
    }
    print OUT "\n";
        
    for(my $y=0;$y<$numYHeaders;$y++) {
        my $yHeader=$inc2header->{ y }->{$y};
        my $yHeaderObject=getHeaderObject($yHeader);
        my $yHeaderMidpoint=$yHeaderObject->{ midpoint };
        
        my %tmpMatrixRow=();
        
        for(my $x=0;$x<$numXHeaders;$x++) {
            my $xHeader=$inc2header->{ x }->{$x};    
            my $xHeaderObject=getHeaderObject($xHeader);
            my $xHeaderMidpoint=$xHeaderObject->{ midpoint };
            
            my $interactionDistance=getInteractionDistance($matrixObject,$yHeaderObject,$xHeaderObject,1);
            next if($interactionDistance == -1);
            next if($interactionDistance > $maxDistance);
            $interactionDistance *= -1 if($xHeaderMidpoint < $yHeaderMidpoint);
            
            my $cScore=$missingValue;
            $cScore = $matrix->{$y}->{$x} if(exists($matrix->{$y}->{$x}));
            $cScore = "NA" if(($cScore eq "") or ($cScore =~ /^NULL$/i) or ($cScore =~ /^NA$/i) or ($cScore =~ /inf$/i) or ($cScore =~ /^nan$/i) or ($cScore eq "NA"));
            $cScore = sprintf "%.1f", $cScore if($cScore ne "NA");
            
            my $binOffset="NA";
            if($matrixObject->{ symmetrical }) {
                $binOffset=int(ceil((($interactionDistance-($headerSizing-$headerSpacing))/$headerSpacing)));            
            } else {
                $binOffset=($x-$y);
            }
            
            next if($binOffset eq "NA");
            
            $tmpMatrixRow{$binOffset}=$cScore;
        }
        
        print OUT "$yHeader";
        for(my $i=-$nOffsetBins;$i<=$nOffsetBins;$i++) {
            my $cScore="NA";
            $cScore = $tmpMatrixRow{$i} if(defined($tmpMatrixRow{$i}));
            print OUT "\t$cScore";
        }
        print OUT "\n";
        
    }
}

my %options;
my $results = GetOptions( \%options,'inputMatrix|i=s','verbose|v','output|o=s','minDistance|minDist=i','maxDistance|maxDist=i') or croak help();
my ($ret,$inputMatrix,$verbose,$output,$minDistance,$maxDistance)=check_options( \%options );

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

#read Matrix
print STDERR "subsetting matrix\n" if($verbose);
my ($subsetMatrix,$subsetSpan)=subsetData($inputMatrix,$header2inc,$minDistance,$maxDistance);
print STDERR "\t$subsetSpan\n" if($verbose);

print STDERR "\n" if($verbose);

my $subsetMatrixFile=$output.".subset.matrix.gz";
print STDERR "Writing matrix to file ($subsetMatrixFile)...\n" if($verbose);
writeMatrix($subsetMatrix,$inc2header,$subsetMatrixFile,$missingValue);
print STDERR "\tcomplete\n" if($verbose);

print STDERR "\n" if($verbose);

my $subsetMatrixObject=getMatrixObject($subsetMatrixFile,$output,$verbose);

$subsetSpan=50;
print STDERR "calculating stacked matrix...\n" if($verbose);
matrix2stacked($subsetMatrixObject,$subsetMatrix,$minDistance,$maxDistance,$missingValue,$subsetSpan);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);