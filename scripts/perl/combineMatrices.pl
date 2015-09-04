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

    my ($inputMatrixArray,$verbose,$output,$combineMode,$tmpDir);
    
    if( exists($opts->{ inputMatrixArray }) ) {
        $inputMatrixArray = $opts->{ inputMatrixArray };
    } else {
        print STDERR "\nERROR: Option inputMatrixArray|i is required.\n";
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
        print STDERR "\nERROR: Option output|o is required.\n";
        help();
    }
    
    if( exists($opts->{ combineMode }) ) {
        $combineMode = $opts->{ combineMode };
    } else {
        $combineMode = "sum";
    }
    
    if( exists($opts->{ tmpDir }) ) {
        $tmpDir = $opts->{ tmpDir };
        $tmpDir =~ s/\/$//;
    } else {
        $tmpDir = "/tmp";
    }
    
    return($inputMatrixArray,$verbose,$output,$combineMode,$tmpDir);
}

sub intro() {
    print STDERR "\n";
    
    print STDERR "Tool:\t\tcombineMatrices.pl\n";
    print STDERR "Version:\t".$cworld::dekker::VERSION."\n";
    print STDERR "Summary:\tcombine matrices [sum,mean,median,min,max]\n";
    
    print STDERR "\n";
}

sub help() {
    intro();
    
    print STDERR "Usage: perl combineMatrices.pl [OPTIONS] -i <inputMatrix>\n";
    
    print STDERR "\n";
    
    print STDERR "Required:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-i@", "[]", "input matrix file array [MULTIPLE]");
    
    print STDERR "\n";
    
    print STDERR "Options:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-v", "[]", "FLAG, verbose mode");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-o", "[]", "prefix for output file(s)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--n", "[]", "optional, prefix for output file");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--cm", "[]", "optional, combine mode [sum,mean,median,min,max]");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--tmp", "[/tmp]", "optional, tmp direction for tmp files");

    print STDERR "\n";
    
    print STDERR "Notes:";
    print STDERR "
    This script can combine multiple matrices.
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

sub combineMatrices($$$$;$) {
    # required
    my $matrixString=shift;
    my $output=shift;
    my $combineMode=shift;
    my $tmpDir=shift;
    # optional
    my $verbose=0;
    $verbose=shift if @_;
    
    my @matrixFileArray=split(/,/,$matrixString);
    
    my %combinedMatrix=();
    
    my @pairwiseFileHandleArray=();
    
    my $numLines=-1;
    my $lastNumLines=-1;
    
    $tmpDir=createTmpDir($tmpDir);
    
    for(my $i=0;$i<@matrixFileArray;$i++) {
        my $inputMatrix=$matrixFileArray[$i];
        
        # get matrix information
        my $matrixObject=getMatrixObject($inputMatrix);
        my $inc2header=$matrixObject->{ inc2header };
        my $header2inc=$matrixObject->{ header2inc };
        my $numYHeaders=$matrixObject->{ numYHeaders };
        my $numXHeaders=$matrixObject->{ numXHeaders };
        my $missingValue=$matrixObject->{ missingValue };
        my $output=$matrixObject->{ output };
        
        my $excludeCis=0;
        my $excludeTrans=0;
        
        print STDERR "\t$output\n" if($verbose);
        
        my $pairwiseFile=matrix2pairwise($matrixObject,$inputMatrix,$excludeCis,$excludeTrans);
        $numLines=getNumberOfLines($pairwiseFile);
        my $pairwiseFileName=getFileName($pairwiseFile);
        
        print STDERR "\t$numLines interactions\n" if($verbose);
        
        croak "pairwise files are not equal length!" if(($lastNumLines != -1) and ($lastNumLines != $numLines));
        
        print STDERR "\tsorting ... " if($verbose);
        my $sortedPairwiseFile=$tmpDir.$pairwiseFileName.".sorted.pairwise.txt.gz";
        system("gunzip -c '".$pairwiseFile."' | grep -v '^#' | sort -k1,1 -k2,2 | gzip > '".$sortedPairwiseFile."'");
        print STDERR "done\n" if($verbose);
        
        system("rm '$pairwiseFile'");
        
        print STDERR "\n" if($verbose);
        
        open(my $tmpFileHandle,inputWrapper($sortedPairwiseFile)) or croak "Could not open file [$sortedPairwiseFile] - $!";
        push(@pairwiseFileHandleArray,$tmpFileHandle);
        
        $lastNumLines=$numLines;
    }
    
    #since all matrices are same size, use info from first
    my $matrixObject=getMatrixObject($matrixFileArray[0]);
    my $inc2header=$matrixObject->{ inc2header };
    my $header2inc=$matrixObject->{ header2inc };
    my $numYHeaders=$matrixObject->{ numYHeaders };
    my $numXHeaders=$matrixObject->{ numXHeaders };

    my $combinedMatrixSum=0;
    
    my $progressBucketSize=ceil($numLines / 1000);
    
    my $pcComplete=0;
    my $lineNum=0;
    my $continue = 1;
    while($continue) {
        
        my @tmpScores=();
        my $lastY="NA";
        my $lastX="NA";
        my $NAFlag=0;
        
        for(my $i=0;$i<@pairwiseFileHandleArray;$i++) {
            my $tmpFileHandle=$pairwiseFileHandleArray[$i];
            my $tmpLine=<$tmpFileHandle>;
            chomp($tmpLine);
            next if(($tmpLine eq "") or ($tmpLine =~ m/^#/));
            
            my @tmp=split(/\t/,$tmpLine);
            my $yHeader=$tmp[0];
            my $y=$header2inc->{ y }->{$yHeader};
            my $xHeader=$tmp[1];
            my $x=$header2inc->{ x }->{$xHeader};
            my $cScore=$tmp[2];
            
            push(@tmpScores,$cScore) if($cScore ne "NA");
            $NAFlag=1 if($cScore eq "NA");
            
            $continue = 0 if(eof($tmpFileHandle));
            next if(eof($tmpFileHandle));
            
            croak "file out of order! [Y - $y - $lastY] [$tmpLine]" if(($lastY ne "NA") and ($lastY != $y));
            croak "file out of order! [X - $x - $lastX] [$tmpLine]" if(($lastX ne "NA") and ($lastX != $x));
            
            $lastY=$y;
            $lastX=$x;
        }
        
        next if(($lastY eq "NA") or ($lastX eq "NA"));
        
        my $tmpScoreStats=listStats(\@tmpScores) if(@tmpScores > 0);
        
        my $combinedScore="NA";      
        croak "combineMode [$combineMode] is invalid!" if((!exists($tmpScoreStats->{$combineMode})) and (@tmpScores > 0));
        $combinedScore=$tmpScoreStats->{$combineMode} if(exists($tmpScoreStats->{$combineMode}));
        
        $combinedScore="NA" if(($NAFlag) and ($combineMode eq "sum"));
        print STDERR "\tWARNING - summing across non equal numbers per interaction!\n" if((($NAFlag) and ($combineMode eq "sum")) and ($combinedScore ne "NA"));
        
        $combinedMatrix{$lastY}{$lastX}=$combinedScore if($combinedScore ne "NA");
        
        $combinedMatrixSum += $combinedScore if($combinedScore ne "NA");
        
        if(($lineNum % $progressBucketSize) == 0) {
            $pcComplete = 100 if($lineNum >= $numLines);
            print STDERR "\e[A" if(($verbose) and ($lineNum != 0));
            printf STDERR "\t%.2f%% complete ($lineNum/".$numLines.")...\n", $pcComplete if($verbose);
            $pcComplete = round((($lineNum/$numLines)*100),2);
        }
        $lineNum++;
    }
    
    for(my $i=0;$i<@pairwiseFileHandleArray;$i++) {
        my $tmpFileHandle=$pairwiseFileHandleArray[$i];
        close($tmpFileHandle);
    }
    
    removeTmpDir($tmpDir);
    
    return(\%combinedMatrix,$combinedMatrixSum);
    
}    
    
my %options;
my $results = GetOptions( \%options,'inputMatrixArray|i=s@','verbose|v','output|o=s','combineMode|cm=s','tmpDir|tmp=s') or croak help();

my ($inputMatrixArray,$verbose,$output,$combineMode,$tmpDir)=check_options( \%options );

intro() if($verbose);

#get the absolute path of the script being used
my $cwd = getcwd();
my $fullScriptPath=abs_path($0);
my @fullScriptPathArr=split(/\//,$fullScriptPath);
@fullScriptPathArr=@fullScriptPathArr[0..@fullScriptPathArr-3];
my $scriptPath=join("/",@fullScriptPathArr);

print STDERR "inputMatrixArray (-i)\t".@{$inputMatrixArray}." files...\n" if($verbose);

print STDERR "\n" if($verbose);

# validate all files are the same
print STDERR "validating identical matrices...\n" if($verbose);
for(my $i1=0;$i1<@{$inputMatrixArray};$i1++) {
    my $inputMatrix_1 = $inputMatrixArray->[$i1];
    
    croak "inputMatrix [$inputMatrix_1] does not exist" if(!(-e $inputMatrix_1));
    
    for(my $i2=0;$i2<@{$inputMatrixArray};$i2++) {
        my $inputMatrix_2 = $inputMatrixArray->[$i2];
        croak "inputMatrix [$inputMatrix_2] does not exist" if(!(-e $inputMatrix_2));
        
        validateIdenticalMatrixStructure($inputMatrix_1,$inputMatrix_2);
    }
}
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

my $matrixString="";
for(my $i=0;$i<@{$inputMatrixArray};$i++) {
    my $inputMatrix = $inputMatrixArray->[$i];
   
    $matrixString .= ",".$inputMatrix if($matrixString ne "");
    $matrixString = $inputMatrix if($matrixString eq "");
    
}

print STDERR "combining matrices...\n" if($verbose);
my ($combinedMatrix,$combinedMatrixSum)=combineMatrices($matrixString,$output,$combineMode,$tmpDir,$verbose);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

# get matrix information
my $matrixObject=getMatrixObject($inputMatrixArray->[0],$output,$verbose);
my $inc2header=$matrixObject->{ inc2header };
my $header2inc=$matrixObject->{ header2inc };
my $numYHeaders=$matrixObject->{ numYHeaders };
my $numXHeaders=$matrixObject->{ numXHeaders };
my $missingValue=$matrixObject->{ missingValue };
my $symmetrical=$matrixObject->{ symmetrical };
my $inputMatrixName=$matrixObject->{ inputMatrixName };
$output=$matrixObject->{ output };

my $combinedMatrixFile=$output.".".$combineMode.".matrix.gz";
print STDERR "Writing matrix to file ($combinedMatrixFile)...\n" if($verbose);
writeMatrix($combinedMatrix,$inc2header,$combinedMatrixFile,"NA");
print STDERR "\tcomplete\n" if($verbose);

print STDERR "\n" if($verbose);

print STDERR "$combinedMatrixFile\t$combinedMatrixSum\n" if($verbose);

print STDERR "\n" if($verbose);