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
    
    my ($inputMatrix,$verbose,$output,$extractBy,$extractCisOnly,$optionalZoomCoordinateArray,$tmpDir);
 
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
    
    if( exists($opts->{ extractBy }) ) {
        $extractBy = $opts->{ extractBy };
        croak "incorrect extractBy [$extractBy] value (chr,region,group)" if(($extractBy ne "chr") and ($extractBy ne "region") and ($extractBy ne "group"));
    } else {
        $extractBy="chr";
    }
    
    if( exists($opts->{ extractCisOnly }) ) {
        $extractCisOnly = 1;
    } else {
        $extractCisOnly = 0;
    }
    
    if( exists($opts->{ optionalZoomCoordinateArray }) ) {
        $optionalZoomCoordinateArray = $opts->{ optionalZoomCoordinateArray };
    } else {
        $optionalZoomCoordinateArray=[];
    }
    
    if( exists($opts->{ tmpDir }) ) {
        $tmpDir = $opts->{ tmpDir };
        $tmpDir =~ s/\/$//;
    } else {
        $tmpDir = "/tmp";
    }
    
    $ret->{ inputMatrix }=$inputMatrix;
    $ret->{ verbose }=$verbose;
    $ret->{ output }=$output;
    $ret->{ extractBy }=$extractBy;
    $ret->{ extractCisOnly }=$extractCisOnly;
    $ret->{ optionalZoomCoordinateArray }=$optionalZoomCoordinateArray;
    $ret->{ tmpDir }=$tmpDir;
    
    return($ret,$inputMatrix,$verbose,$output,$extractBy,$extractCisOnly,$optionalZoomCoordinateArray,$tmpDir);
}

sub getYHeaders($$;$) {
    # required
    my $datasetFile=shift;
    my $extractBy=shift;
    # optional
    my $verbose=0;
    $verbose=shift if @_;
    
    
    my $headers={};
    my $header2index={};
    my $headerSizes={};
    my $headerBoundaries={};
    
    my $line="";
    my $lineNum=0;
    my $lastSubMatrix="";
    
    my $xHeaderLine="";
    
    my @rowBlockFiles=();
    
    my $nLines = getNumberOfLines($datasetFile);
    my $progressBucketSize=ceil($nLines / 1000);
    my $pcComplete=0;
    
    open(IN,inputWrapper($datasetFile)) or croak "Could not open file [$datasetFile] - $!";
    while($line = <IN>) {
        chomp($line);
        next if(($line eq "") or ($line =~ m/^#/));
        
        if($lineNum == 0) {
            $xHeaderLine=$line;
        } else {
            my $head=(split(/\t/,$line))[0];
            my $ySubMatrix=header2subMatrix($head,$extractBy);
            
            my $yIndex=0;
            $yIndex=$headerSizes->{ y }->{$ySubMatrix} if(exists($headerSizes->{ y }->{$ySubMatrix}));
            
            $headers->{ y }->{$ySubMatrix}->{$yIndex}=$head;
            $header2index->{ y }->{$ySubMatrix}->{$head}=$yIndex;
            
            $headerSizes->{ y }->{$ySubMatrix}++;
        } 
        $pcComplete = 100 if($lineNum == ($nLines-1));
        print STDERR "\e[A" if(($verbose) and ($lineNum != 0));
        printf STDERR "\t%.2f%% complete ($lineNum/$nLines)...\n", $pcComplete if($verbose);
        $pcComplete = round((($lineNum/$nLines)*100),2);
        $lineNum++;
    }
    close(IN);
    
    $pcComplete=100;
    print STDERR "\e[A" if($verbose);
    printf STDERR "\t%.2f%% complete ($lineNum/$nLines)...\n", $pcComplete if($verbose);
    
    return($headers,$header2index,$headerSizes);
}

sub extractRowBlocks($$$$$$$) {
    my $datasetFile=shift;
    my $datasetFileName=shift;
    my $extractBy=shift;
    my $blockMatrix=shift;
    my $extractCisOnly=shift;
    my $zoomData=shift;
    my $goodZoomFlag=shift;
    
    my $line="";
    my $lineNum=0;
    my $lastSubMatrix="";
    my $suffix="";
    my $file="";
    my $xHeaderLine="";
    
    my $rowBlockFiles={};
    
    my $suffixInit={};
    
    open(IN,inputWrapper($datasetFile)) or croak "Could not open file [$datasetFile] - $!";
    while($line = <IN>) {
        chomp($line);
        next if(($line eq "") or ($line =~ m/^#/));
        
        if($lineNum == 0) {
            my $matrixSize=(split(/\t/,$line))[0];
            $xHeaderLine=$line;
            $xHeaderLine =~ s/$matrixSize//;
        } else {
            
            my $header=(split(/\t/,$line))[0];
            
            my $subMatrix="";
            $subMatrix=header2subMatrix($header,$extractBy);

            my $optionalZoomCoordinateName="";
            my $headerObject=getHeaderObject($header);
            my $headerChromosome=$headerObject->{ chromosome };
            my $headerStart=$headerObject->{ start };
            my $headerEnd=$headerObject->{ end };
            
            my $keepLine=0;
            if($goodZoomFlag) {
                foreach my $z ( keys %{$zoomData} ) {                                                
                    next if($keepLine);
                    my $tmpZoomData=$zoomData->{$z};
                    $keepLine = 1 if(($headerChromosome eq $tmpZoomData->{ chromosome }) and (isOverlapping($tmpZoomData->{ start },$tmpZoomData->{ end },$headerStart,$headerEnd)));
                    $optionalZoomCoordinateName="__".$tmpZoomData->{ chromosome }."_".$tmpZoomData->{ start }."_".$tmpZoomData->{ end } if($keepLine);
                }
            }
            
            $subMatrix .= $optionalZoomCoordinateName if($extractCisOnly);
            $subMatrix = "userzoom" if(($goodZoomFlag) and ($extractCisOnly == 0));
            
            croak "bad subMatrix ($subMatrix) - line# $lineNum" if($subMatrix eq "");
            
            my $interactionType="all";
            $interactionType="cis" if(($subMatrix eq $blockMatrix) and ($blockMatrix ne "all"));
            $interactionType="trans" if(($subMatrix ne $blockMatrix) and ($blockMatrix ne "all"));
            
            # only proces CIS if -eco enabled
            next if(($interactionType eq "trans") and ($extractCisOnly));
            
            next if(($goodZoomFlag) and ($keepLine == 0));
            
            # if selected - process cis only
            $lastSubMatrix=$subMatrix if(($interactionType eq "trans") and ($extractCisOnly));
            $lineNum++ if(($interactionType eq "trans") and ($extractCisOnly));
                            
            if($subMatrix ne $lastSubMatrix) { #close and open new file
                if($lastSubMatrix ne "") {
                    close(TMP);
                }
                
                $suffix=$subMatrix."__".$interactionType."__all" if($blockMatrix eq "");
                $suffix=$subMatrix."__".$blockMatrix."__".$interactionType if(($blockMatrix ne "") and ($subMatrix ne $blockMatrix));
                $suffix=$subMatrix."__".$interactionType if(($blockMatrix ne "") and ($subMatrix eq $blockMatrix));
                
                $file=$datasetFileName."___".$suffix.".matrix.gz";
                
                croak "bad filename [$suffix] - blockMatrix [$blockMatrix]" if($suffix eq "");
                
                open(TMP,outputWrapper($file,"",1)) or croak "Could not open file [$file] - $!";
                $rowBlockFiles->{$suffix}->{ file }=$file;
                $rowBlockFiles->{$suffix}->{ subMatrix}=$subMatrix;
                print TMP "$xHeaderLine\n" if(!exists($suffixInit->{$suffix}));
                $suffixInit->{$suffix}=1;
            }
            print TMP "$line\n";
            
            $lastSubMatrix=$subMatrix;
        }
        $lineNum++;
    }
    close(IN);
    close(TMP);
    
    return($rowBlockFiles);
}


sub transposeFiles($$;$) {
    # required
    my $fileList=shift;
    my $matrixName=shift;
    # optional
    my $verbose=0;
    $verbose = shift if @_;
    
    my $transposedFiles={};
    
    foreach my $suffix ( keys %$fileList ) {    
        my $file=$fileList->{$suffix}->{ file };
        my $subMatrix=$fileList->{$suffix}->{ subMatrix };
        
        # now transpose the matrix
        print STDERR "\ttransposing [$file] ...\n" if($verbose);
        
        my $transposedMatrix=transposeMatrix($file,$matrixName."__".getSmallUniqueString().".transpose");
        
        $transposedFiles->{$suffix}->{ file } = $transposedMatrix;
        $transposedFiles->{$suffix}->{ subMatrix } = $subMatrix;
        
    }
    
    return($transposedFiles);
}

sub extractSubMatrices($$$$$$$) {
    my $fileList=shift;
    my $dataFileName=shift;
    my $output=shift;
    my $extractBy=shift;
    my $extractCisOnly=shift;
    my $zoomData=shift;
    my $goodZoomFlag=shift;

    my @files=();
    
    foreach my $suffix ( keys %$fileList ) {    
        my $file=$fileList->{$suffix}->{ file };    
        my $blockMatrix=$fileList->{$suffix}->{ subMatrix };
        
        my $rowBlockFiles=extractRowBlocks($file,$dataFileName,$extractBy,$blockMatrix,$extractCisOnly,$zoomData,$goodZoomFlag);
          
        foreach my $suffix ( keys %$rowBlockFiles ) {    
            
            my $tmp_matrix=$rowBlockFiles->{$suffix}->{ file };
            my $subMatrix=$rowBlockFiles->{$suffix}->{ subMatrix };
            my $tmp_transposedMatrix=transposeMatrix($tmp_matrix,$output."___".$suffix);
        }
        
    }
    
}    

sub removeFiles($) {
    my $fileList=shift;
    
    foreach my $suffix ( keys %$fileList ) {    
        my $file=$fileList->{$suffix}->{ file };
        my $subMatrix=$fileList->{$suffix}->{ subMatrix };
        
        system("rm '".$file."'");
        
    }
}    

sub intro() {
    print STDERR "\n";
    
    print STDERR "Tool:\t\t".$tool."\n";
    print STDERR "Version:\t".$cworld::dekker::VERSION."\n";
    print STDERR "Summary:\tthis script can extract sub matrices classified by chr/group/name\n";
    
    print STDERR "\n";
}

sub help() {
    intro();
    
    print STDERR "Usage: perl extractSubMatrices.pl [OPTIONS] -i <inputMatrix>\n";
    
    print STDERR "\n";
    
    print STDERR "Required:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-i", "[]", "input matrix file");
    
    print STDERR "\n";
    print STDERR "Options:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-v", "[]", "FLAG, verbose mode");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-o", "[]", "prefix for output file(s)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ozc", "[]", "@ [MULTIPLE], genomic coordinates for subset (-ozc chr:start-end)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--eb", "[chr]", "optional, classifier by which to extract sub matrices (chr,region,group)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--eco", "[]", "FLAG, extract only cis matrices");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--usn", "[]", "FLAG, use short input file names");    
    printf STDERR ("\t%-10s %-10s %-10s\n", "--tmp", "[/tmp]", "optional, tmp direction for tmp files");
    print STDERR "\n";
    
    print STDERR "Notes:";
    print STDERR "
    This script extracts all sub matrices of a given input matrix.
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

my %options;
my $results = GetOptions( \%options,'inputMatrix|i=s','verbose|v','output|o=s','extractBy|eb=s','extractCisOnly|eco','optionalZoomCoordinateArray|ozc=s@','tmpDir|tmp=s') or croak help();
my ($ret,$inputMatrix,$verbose,$output,$extractBy,$extractCisOnly,$optionalZoomCoordinateArray,$tmpDir)=check_options( \%options );

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

$tmpDir=createTmpDir($tmpDir);

# set up tmpDir prior to optional job dir / name etc
my $tmpDirInputMatrixName=$tmpDir.$inputMatrixName;

print STDERR "\n" if($verbose);

# initial pass on matrix
print STDERR "parsing initial matrix...\n" if($verbose);
my ($headers,$header2index,$headerSizes)=getYHeaders($inputMatrix,$extractBy); 
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

my $zoomData={};
my $goodZoomFlag=0;
for(my $z=0;$z<@{$optionalZoomCoordinateArray};$z++) {
    my $optionalZoomCoordinate=$optionalZoomCoordinateArray->[$z];
    next if(!(validateZoomCoordinate($optionalZoomCoordinate)));
    
    my ($tmpZoomData)=splitCoordinate($optionalZoomCoordinate);
    $zoomData->{$z}=$tmpZoomData;
    print STDERR "optionalZoomCoordinates\t".$tmpZoomData->{ chromosome }.":".$tmpZoomData->{ start }."-".$tmpZoomData->{ end }."\n" if($verbose);
    
    $goodZoomFlag = 1 if($tmpZoomData->{ flag });
}

if($goodZoomFlag) {
    print STDERR "\nchanging extract by to (chr)\n" if($verbose);
    $extractBy="chr";
}

# ensure zoom coordinates do not overlap (must be distinct)
for(my $z1=0;$z1<(keys %{$zoomData});$z1++) {
    my $tmpZoomData1=$zoomData->{$z1};
    my $tmpZoomData1Name=$tmpZoomData1->{ name };
    for(my $z2=0;$z2<(keys %{$zoomData});$z2++) {
        my $tmpZoomData2=$zoomData->{$z2};
        my $tmpZoomData2Name=$tmpZoomData2->{ name };
        
        next if($z1 == $z2);
        
        if(isOverlapping($tmpZoomData1->{ start },$tmpZoomData1->{ end },$tmpZoomData2->{ start },$tmpZoomData2->{ end },$tmpZoomData1->{ chromosome },$tmpZoomData2->{ chromosome })) {
            croak "zoom coordinates cannot overlap - must be distinct [$tmpZoomData1Name] [$tmpZoomData2Name]";
        }
    }
}

print STDERR "extracting row blocks ...\n" if($verbose);
my $rowBlockFiles=extractRowBlocks($inputMatrix,$tmpDirInputMatrixName,$extractBy,"all",$extractCisOnly,$zoomData,$goodZoomFlag);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

# transpose all row block files
print STDERR "transposing all row blocks ...\n" if($verbose);
my $transposedFiles=transposeFiles($rowBlockFiles,$tmpDirInputMatrixName,$verbose);

print STDERR "\n" if($verbose);

print STDERR "extracting sub matrices ...\n" if($verbose);
# now re-extract row blocks to get sub matrices
extractSubMatrices($transposedFiles,$tmpDirInputMatrixName,$output,$extractBy,$extractCisOnly,$zoomData,$goodZoomFlag);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

removeFiles($rowBlockFiles);
removeFiles($transposedFiles);
removeTmpDir($tmpDir);
