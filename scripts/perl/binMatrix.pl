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

    my ($inputMatrix,$verbose,$output,$binSize,$binStep,$binMode,$excludeZero,$niceCoordinatesFlag);
    
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
    
    if( exists($opts->{ binSize }) ) {
        $binSize = $opts->{ binSize };
    } else {
        print STDERR "\nERROR: Option binSize|bsize is required.\n";
        help();
    }
    
    if( exists($opts->{ binStep }) ) {
        $binStep = $opts->{ binStep };
    } else {
        print STDERR "\nERROR: Option binStep|bstep is required.\n";
        help();
    }
    
    if( exists($opts->{ binMode }) ) {
        $binMode = $opts->{ binMode };
    } else {
        print STDERR "\nERROR: Option binMode|bmode is required.\n";
        help();
    }
    
    if( exists($opts->{ excludeZero }) ) {
        $excludeZero = 1;
    } else {
        $excludeZero = 0;
    }
    
    if( exists($opts->{ niceCoordinates }) ) {
        $niceCoordinatesFlag = $opts->{ niceCoordinates };
    } else {
        $niceCoordinatesFlag=0;
    }

    
    return($inputMatrix,$verbose,$output,$binSize,$binStep,$binMode,$excludeZero,$niceCoordinatesFlag);
}

sub intro() {
    print STDERR "\n";
    
    print STDERR "Tool:\t\tbinMatrix.pl\n";
    print STDERR "Version:\t".$cworld::dekker::VERSION."\n";
    print STDERR "Summary:\tbin/aggregrate matrix into fixed sized intervals\n";
    
    print STDERR "\n";
}

sub help() {
    intro();
    
    print STDERR "Usage: perl binMatrix.pl [OPTIONS] -i <inputMatrix> -bsize <binSize> -bstep <binStep> -bmode <binMode>\n";
    
    print STDERR "\n";
    
    print STDERR "Required:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-i", "[]", "input matrix file");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-o", "[]", "prefix for output file(s)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--bsize", "[]", "bin size in bp");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--bstep", "[]", "bin step in factor of overlap, 1=non overlapping, 2-inf=overlapping");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--bmode", "[]", "bin mode, how to aggregrate data, mean,sum,median etc.");
    
    print STDERR "\n";
    
    print STDERR "Options:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-v", "[]", "FLAG, verbose mode");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ez", "[]", "FLAG, ignore 0s in all calculations");    
    printf STDERR ("\t%-10s %-10s %-10s\n", "--nc", "[]", "FLAG, use nice(clean) coordinates.  start first bin at position cleanly divisible by window size parameter.");    
    
    print STDERR "\n";
    
    print STDERR "Notes:";
    print STDERR "
    This script can bin a matrix into fixed size intervals and aggregrate the data.
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

sub writeBinPositions($$$$$) {    
    my $matrixObject=shift;
    my $output=shift;
    my $binSize=shift;
    my $binStep=shift;
    my $niceCoordinatesFlag=shift;
    
    my $contigList=$matrixObject->{ contigList };
    my $index2contig=$matrixObject->{ index2contig };
    my $numContigs=$matrixObject->{ numContigs };
    
    my $binOffsets={};
    my $binOffset=0;
    for(my $c=0;$c<$numContigs;$c++) {
        my $contig=$index2contig->{ xy }->{$c};
        my $contigLength=$contigList->{$contig}->{ contigLength };
        $binOffsets->{$contig}=$binOffset;
        $binOffset +=ceil($contigLength/($binSize/$binStep));
    }
    
    my $totalNumBins=0;
    my $maxBins=10000;
    
    my $binHeaderFile=$output.".headers";
    open(OUT,outputWrapper($binHeaderFile)) or croak "Could not open file [$binHeaderFile] - $!";
    for(my $c=0;$c<$numContigs;$c++) {
        my $contig=$index2contig->{ xy }->{$c};
        my $contigLength=$contigList->{$contig}->{ contigLength };
        my $contigStart=$contigList->{$contig}->{ contigStart };
        my $contigEnd=$contigList->{$contig}->{ contigEnd };
        my $contigAssembly=$contigList->{$contig}->{ contigAssembly };
        my $contigChromosome=$contigList->{$contig}->{ contigChromosome };
        my $numBins = ceil($contigLength/($binSize/$binStep));
        my $binOffset=$binOffsets->{$contig};
        
        for(my $i=0;$i<$numBins;$i++) {
            my $binStart = $contigStart+($i*($binSize/$binStep));
            my $binEnd = $binStart + $binSize;
            $binEnd=$contigEnd if($binEnd > $contigEnd);
            my $binIndex=$contig."__".$i;
            my $binName=$binIndex."|".$contigAssembly."|".$contigChromosome.":".$binStart."-".$binEnd;
            print OUT "$binIndex\t$contig\t$contigChromosome\t$binStart\t$binEnd\t$binName\n";
            
            $totalNumBins++;
            croak "attempted binned matrix is too large! ($totalNumBins > $maxBins)" if($totalNumBins > $maxBins);
        }
    }
    close(OUT);
    
    return($binHeaderFile,$binOffsets);
}

sub reOrient($$$$$$) {
    my $yHeader=shift;
    my $yBinIndexStart=shift;
    my $yBinIndexEnd=shift;
    my $xHeader=shift;
    my $xBinIndexStart=shift;
    my $xBinIndexEnd=shift;
    
    if($yBinIndexStart <= $xBinIndexStart) {
        return($yHeader,$yBinIndexStart,$yBinIndexEnd,$xHeader,$xBinIndexStart,$xBinIndexEnd);
    } else {
        return($xHeader,$xBinIndexStart,$xBinIndexEnd,$yHeader,$yBinIndexStart,$yBinIndexEnd);
    }
    
}

sub matrix2binned($$$$$$$$) {
    my $matrixObject=shift;
    my $inputMatrix=shift;
    my $output=shift;
    my $binOffsets=shift;
    my $binSize=shift;
    my $binStep=shift;
    my $binMode=shift;
    my $excludeZero=shift;
    
    my $contigList=$matrixObject->{ contigList };
    my $header2contig=$matrixObject->{ header2contig };
    my $index2contig=$matrixObject->{ index2contig };
    my $numContigs=$matrixObject->{ numContigs };
    my $inc2header=$matrixObject->{ inc2header };
    my $header2inc=$matrixObject->{ header2inc };
    my $numYHeaders=$matrixObject->{ numYHeaders };
    my $numXHeaders=$matrixObject->{ numXHeaders };
    my $missingValue=$matrixObject->{ missingValue };
    my $verbose=$matrixObject->{ verbose };
    
    my $lineNum=0;
    my @xHeaders=();
    
    my $pcComplete=0;
    my $nLines = getNumberOfLines($inputMatrix)-1;
    
    my $binnedPairwiseFile=$output.".bin";
    open(OUT,outputWrapper($binnedPairwiseFile)) or croak "Could not open file [$binnedPairwiseFile] - $!";
    
    my $initFlag=1;
    open(IN,inputWrapper($inputMatrix)) or croak "Could not open file [$inputMatrix] - $!";
    while(my $line = <IN>) {
        chomp($line);
        next if(($line eq "") or ($line =~ m/^#/));
        
        if($lineNum == 0) {
            @xHeaders=split(/\t/,$line);
        } else {
            my @data=split(/\t/,$line);
            my $dsize=@data;
            
            my $yHeader=$data[0];
            
            for(my $d=1;$d<$dsize;$d++) {
            
                my $yContig=$header2contig->{ y }->{$yHeader};
                my $yContigStart=$contigList->{$yContig}->{ contigStart };
                my $yBinOffset=$binOffsets->{$yContig};
                my $yHeaderObject=getHeaderObject($yHeader);
                my $yHeaderMidpoint=$yHeaderObject->{ midpoint }-$yContigStart;
                my $yBinIndexEnd=floor($yHeaderMidpoint/($binSize/$binStep))+$yBinOffset;    
                my $yBinIndexStart = ($yBinIndexEnd-($binStep-1));
                $yBinIndexStart = $yBinOffset if($yBinIndexStart < $yBinOffset);
                
                my $xHeader=$xHeaders[$d];
                
                my $cScore=$data[$d];
                $cScore = "NA" if($cScore !~ (/^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/));
                $cScore = "NA" if(($cScore eq "") or ($cScore =~ /^NULL$/i) or ($cScore =~ /^NA$/i) or ($cScore =~ /inf$/i) or ($cScore =~ /^nan$/i) or ($cScore eq "NA"));
                
                next if(($cScore ne "NA") and ($excludeZero) and ($cScore == 0));
                # if sum - skip 0s
                next if(($cScore ne "NA") and ($binMode eq "sum") and ($cScore == 0));
                
                $cScore = sprintf("%.4f",$cScore) if($cScore ne "NA");
                
                my $xContig=$header2contig->{ x }->{$xHeader};
                my $xContigStart=$contigList->{$xContig}->{ contigStart };
                my $xBinOffset=$binOffsets->{$xContig};
                my $xHeaderObject=getHeaderObject($xHeader);
                my $xHeaderMidpoint=$xHeaderObject->{ midpoint }-$xContigStart;

                my $xBinIndexEnd=floor($xHeaderMidpoint/($binSize/$binStep))+$xBinOffset;    
                my $xBinIndexStart = ($xBinIndexEnd-($binStep-1));
                $xBinIndexStart = $xBinOffset if($xBinIndexStart < $xBinOffset);
                
                my $tmpYHeader=$yHeader;
                my $tmpXHeader=$xHeader;
                ($tmpYHeader,$yBinIndexStart,$yBinIndexEnd,$tmpXHeader,$xBinIndexStart,$xBinIndexEnd)=reOrient($tmpYHeader,$yBinIndexStart,$yBinIndexEnd,$tmpXHeader,$xBinIndexStart,$xBinIndexEnd);
                
                my $headerKey=$tmpYHeader."__".$tmpXHeader;
                
                print OUT "$yBinIndexStart\t$yBinIndexEnd\t$xBinIndexStart\t$xBinIndexEnd\t$cScore\t$headerKey\n";
                
            }
            
            $pcComplete = 100 if($lineNum == ($nLines-1));
            print STDERR "\e[A" if(($verbose) and ($lineNum > 1));
            printf STDERR "\t%.2f%% complete ($lineNum/".($nLines-1).")...\n", $pcComplete if($verbose);
            $pcComplete = round((($lineNum/$nLines)*100),2);
    
        }
                
        $lineNum++;
    }
    close(IN);
    
    close(OUT);
    
    return($binnedPairwiseFile);
}

sub parseBinHeaders($) {
    my $binHeaderFile=shift;
    
    my $inc2bin={};    
    my $bin2inc={};

    my $maxBins=10000;
    
    my $numBins=0;
    open(IN,inputWrapper($binHeaderFile)) or croak "Could not open file [$binHeaderFile] - $!";
    while(my $line = <IN>) {
        chomp($line);
        next if(($line eq "") or ($line =~ m/^#/));
        
        my ($binIndex,$binContig,$binChromosome,$binStart,$binEnd,$binName)=split(/\t/,$line);

        $inc2bin->{ x }->{$numBins}=$binName;
        $bin2inc->{ x }->{$binName}=$numBins;
        $inc2bin->{ y }->{$numBins}=$binName;
        $bin2inc->{ y }->{$binName}=$numBins;
        $inc2bin->{ xy }->{$numBins}=$binName;
        $bin2inc->{ xy }->{$binName}=$numBins;
        $numBins++;

        croak "attempted binned matrix is too large! (> $maxBins)" if($numBins > $maxBins);

    }
    close(IN);
    
    return($inc2bin,$bin2inc);
        
}

sub binnedPairwise2matrix($$$$$$$) {
    # this method reads in the pairwise file, and the does the score aggreration
    my $matrixObject=shift;
    my $binnedPairwiseFile=shift;
    my $inc2bin=shift;
    my $bin2inc=shift;
    my $binSize=shift;
    my $binStep=shift;
    my $binMode=shift;
    
    my $verbose=$matrixObject->{ verbose };
    
    my $numYHeaders=keys(%{$bin2inc->{ y }});
    my $numXHeaders=keys(%{$bin2inc->{ x }});
    
    croak "matrix not symmetrical ($numYHeaders | $numXHeaders)" if($numYHeaders != $numXHeaders);
    
    my $memoryBinStart=0;
    my $lastMemoryBinStart=0;
    
    my $matrix={};
    my $matrixContainer={};
    
    my $interactionCount=0;
    my $nInteractions=ceil((($numYHeaders*$numXHeaders)/2)+$numYHeaders);
    my $progressBucketSize=ceil($nInteractions / 1000);
    my $pcComplete=0;
    
    open(IN,inputWrapper($binnedPairwiseFile)) or croak "Could not open file [$binnedPairwiseFile] - $!";
    while(my $line = <IN>) {
        chomp($line);
        next if(($line eq "") or ($line =~ m/^#/));
        
        my ($yIndexStart,$yIndexEnd,$xIndexStart,$xIndexEnd,$cScore)=split(/\t/,$line);
    
        next if($cScore eq "NA");
        
        if($yIndexStart != $memoryBinStart) {
            for(my $tmpY=$lastMemoryBinStart;$tmpY<$yIndexStart;$tmpY++) {
                # print STDERR "\tbuilding $tmpY row...\n";
                
                for(my $tmpX=0;$tmpX<$numXHeaders;$tmpX++) {
                
                    # only work above diagonal
                    next if($tmpY > $tmpX);
                    
                    my $tmpArr=[];
                    $tmpArr=$matrixContainer->{$tmpY}->{$tmpX} if(exists($matrixContainer->{$tmpY}->{$tmpX}));
                    my $nScores = scalar @$tmpArr;
                    
                    my $aggregrateScore="NA";
                    
                    if($nScores > 0) {
                        my $tmpArrStats=listStats($tmpArr);
                        $aggregrateScore=$tmpArrStats->{ $binMode };
                    }
                    
                    # now clear out the data from the container
                    delete $matrixContainer->{$tmpY}->{$tmpX};
                    
                    # now log the score in the matrix
                    #print STDERR "\tbuilding $tmpY,$tmpX | $tmpX,$tmpY == $aggregrateScore\n";
                    croak "matrix y,x $tmpY,$tmpX already exists (".$matrix->{$tmpY}->{$tmpX}." ? $aggregrateScore)" if(exists($matrix->{$tmpY}->{$tmpX}));
                    croak "matrix x,y $tmpX,$tmpY already exists (".$matrix->{$tmpX}->{$tmpY}." ? $aggregrateScore)" if(exists($matrix->{$tmpX}->{$tmpY}));
                    
                    $matrix->{$tmpY}->{$tmpX}=$aggregrateScore;
                    $matrix->{$tmpX}->{$tmpY}=$aggregrateScore;
                    
                    if(($interactionCount % $progressBucketSize) == 0) {
                        $pcComplete = 100 if($interactionCount == ($interactionCount-1));
                        print STDERR "\e[A" if(($verbose) and ($interactionCount != 0));
                        printf STDERR "\t%.2f%% complete ($interactionCount/$nInteractions)...\n", $pcComplete if($verbose);
                        $pcComplete = round((($interactionCount/$nInteractions)*100),2);
                    }
                    $interactionCount++;

                }
                
            }
            $lastMemoryBinStart=$yIndexStart;
        }
            
        $memoryBinStart=$yIndexStart;        

        for(my $y=$yIndexStart;$y<=$yIndexEnd;$y++) {
            for(my $x=$xIndexStart;$x<=$xIndexEnd;$x++) {
                
                # only work above diagonal
                next if($y > $x);
                
                # running tally of scores
                push(@{$matrixContainer->{$y}->{$x}},$cScore);
            }
        }
        
    }
    close(IN);
    
    # clear out previous indices
    for(my $tmpY=$lastMemoryBinStart;$tmpY<$numYHeaders;$tmpY++) {
        for(my $tmpX=0;$tmpX<$numXHeaders;$tmpX++) {
        
            # only work above diagonal
            next if($tmpY > $tmpX);
            
            my $tmpArr=[];
            $tmpArr=$matrixContainer->{$tmpY}->{$tmpX} if(exists($matrixContainer->{$tmpY}->{$tmpX}));
            my $nScores = scalar @$tmpArr;
            
            my $aggregrateScore="NA";
            
            if($nScores > 0) {
                my $tmpArrStats=listStats($tmpArr);
                $aggregrateScore=$tmpArrStats->{ $binMode };
            }
            
            # now clear out the data from the container
            delete $matrixContainer->{$tmpY}->{$tmpX};
            
            # now log the score in the matrix
            croak "matrix y,x $tmpY,$tmpX already exists (".$matrix->{$tmpY}->{$tmpX}." ? $aggregrateScore)" if(exists($matrix->{$tmpY}->{$tmpX}));
            croak "matrix x,y $tmpX,$tmpY already exists (".$matrix->{$tmpX}->{$tmpY}." ? $aggregrateScore)" if(exists($matrix->{$tmpX}->{$tmpY}));
            
            $matrix->{$tmpY}->{$tmpX}=$aggregrateScore;
            $matrix->{$tmpX}->{$tmpY}=$aggregrateScore;
        }
        
    }    
    
    $pcComplete=100;
    print STDERR "\e[A" if($verbose);
    printf STDERR "\t%.2f%% complete ($nInteractions/$nInteractions)...\n", $pcComplete if($verbose);
    
    return($matrix);
}

sub makeNiceContigs($$$) {
    my $matrixObject=shift;
    my $binSize=shift;
    my $binStep=shift;

    my $contigList=$matrixObject->{ contigList };
    my $index2contig=$matrixObject->{ index2contig };
    my $numContigs=$matrixObject->{ numContigs };
    
    for(my $c=0;$c<$numContigs;$c++) {
        my $contig=$index2contig->{$c};
        my $contigLength=$contigList->{$contig}->{ contigLength };
        my $contigStart=$contigList->{$contig}->{ contigStart };
        my $contigEnd=$contigList->{$contig}->{ contigEnd };
        my $contigAssembly=$contigList->{$contig}->{ contigAssembly };
        my $contigChromosome=$contigList->{$contig}->{ contigChromosome };
        my $numBins = ceil($contigLength/($binSize/$binStep));        
        
        my $contigStartOffset=($contigStart % $binSize)-1;
        $contigStart-=$contigStartOffset;
        
        my $contigEndOffset=$binSize-($contigEnd % $binSize)+1;
        $contigEnd+=$contigEndOffset;
        
        $contigList->{$contig}->{ contigStart }=$contigStart;
        $contigList->{$contig}->{ contigEnd }=$contigEnd;
        $contigList->{$contig}->{ contigLength }=($contigEnd-$contigStart);
    }
    
    $matrixObject->{ contigList }=$contigList;
    
    return($matrixObject);
}
            
my %options;
my $results = GetOptions( \%options,'inputMatrix|i=s','verbose|v','output|o=s','binSize|bsize=s','binStep|bstep=s','binMode|bmode=s','excludeZero|ez','niceCoordinates|nc') or croak help();

my ($inputMatrix,$verbose,$output,$binSize,$binStep,$binMode,$excludeZero,$niceCoordinatesFlag) = check_options( \%options );

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
my $inputMatrixName=$matrixObject->{ inputMatrixName };
$output=$matrixObject->{ output };

# append out suffix
$output .= "__".$binSize."__".$binStep;

$matrixObject=makeNiceContigs($matrixObject,$binSize,$binStep) if($niceCoordinatesFlag);

print STDERR "writing bin positions ... \n" if($verbose);
my ($binHeaderFile,$binOffsets)=writeBinPositions($matrixObject,$output,$binSize,$binStep,$niceCoordinatesFlag);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

print STDERR "matrix2bin ... \n" if($verbose);
my ($binnedPairwiseFile)=matrix2binned($matrixObject,$inputMatrix,$output,$binOffsets,$binSize,$binStep,$binMode,$excludeZero);
print STDERR "\tdone\n" if($verbose);

system("sort -k1,1n -k2,2n -k3,3n -k4,4n $binnedPairwiseFile -o $binnedPairwiseFile");

print STDERR "\n" if($verbose);

print STDERR "parsing headers ... \n" if($verbose);
my ($inc2bin,$bin2inc)=parseBinHeaders($binHeaderFile);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

my $numBins=keys %{$inc2bin->{ xy }};
print STDERR "building binned matrix $numBins x $numBins\n" if($verbose);

print STDERR "\n" if($verbose);

print STDERR "building matrix ... \n" if($verbose);
my ($matrix)=binnedPairwise2matrix($matrixObject,$binnedPairwiseFile,$inc2bin,$bin2inc,$binSize,$binStep,$binMode);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

my $binnedMatrixFile=$output.".matrix.gz";
print STDERR "writing matrix ... \n" if($verbose);
writeMatrix($matrix,$inc2bin,$binnedMatrixFile,"NA");
print STDERR "\tdone\n" if($verbose);

system("rm $binnedPairwiseFile");
system("rm $binHeaderFile");

print STDERR "\n" if($verbose);