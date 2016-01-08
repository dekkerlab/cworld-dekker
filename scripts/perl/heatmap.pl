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

use GD::Simple;

my $tool=(split(/\//,abs_path($0)))[-1];

use cworld::dekker;

# def global for PI
use constant PI    => 4 * atan2(1, 1);

sub check_options {
    my $opts = shift;

    my ($inputMatrixArray,$verbose,$output,$outputSuffix,$logTransform,$imageSize,$pixelSize,$y_pixelSize,$x_pixelSize,$maxImageDim,$drawPixelBorder,$omitContigBorder,$drawLabel,$drawScores,$colorScaleStart,$colorScaleEnd,$colorScaleStartTile,$colorScaleEndTile,$posColorString,$negColorString,$missingColor,$highlightColor,$elementBedFile,$imageQuality,$scaleFragmentSizes,$drawTriangle,$drawDiamond,$transparentBGFlag,$embed_meta,$contigSpacing,$scaleMode,$transparency);
    
    my $ret={};
    
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
        $output = "";
    }
    
    if( exists($opts->{ outputSuffix }) ) {
        $outputSuffix = $opts->{ outputSuffix };
    } else {
        $outputSuffix = "";
    }
    
    if( exists($opts->{ logTransform }) ) {
        $logTransform = $opts->{ logTransform };
    } else {
        $logTransform = 0;
    }
    
    if( exists($opts->{ imageSize }) ) {
        $imageSize = $opts->{ imageSize };
    } else {
        $imageSize = 800;
    }
    
    if( exists($opts->{ pixelSize }) ) {
        $pixelSize = $opts->{ pixelSize };
    } else {
        $pixelSize = "NA";
    }
    
    if( exists($opts->{ y_pixelSize }) ) {
        $y_pixelSize = $opts->{ y_pixelSize };
    } else {
        $y_pixelSize = "NA";
    }
    
    if( exists($opts->{ x_pixelSize }) ) {
        $x_pixelSize = $opts->{ x_pixelSize };
    } else {
        $x_pixelSize ="NA";
    }
    
    if( exists($opts->{ maxImageDim }) ) {
        $maxImageDim = $opts->{ maxImageDim };
    } else {
        $maxImageDim = 32000;
    }
    
    if( exists($opts->{ drawPixelBorder }) ) {
        $drawPixelBorder = 1;
    } else {
        $drawPixelBorder=0;
    }
    
    if( exists($opts->{ omitContigBorder }) ) {
        $omitContigBorder = 1;
    } else {
        $omitContigBorder=0;
    }
    
    
    if( exists($opts->{ drawLabel }) ) {
        $drawLabel = 1;
    } else {
        $drawLabel = 0;
    }
    
    if( exists($opts->{ drawScores }) ) {
        $drawScores = 1;
    } else {
        $drawScores = 0;
    }
    
    if( exists($opts->{ colorScaleStart }) ) {
        $colorScaleStart = $opts->{ colorScaleStart };
    } else {
        $colorScaleStart = "NA";
    }
    
    if( exists($opts->{ colorScaleEnd }) ) {
        $colorScaleEnd = $opts->{ colorScaleEnd };
    } else {
        $colorScaleEnd = "NA";
    }
    
    if( exists($opts->{ colorScaleStartTile }) ) {
        $colorScaleStartTile = $opts->{ colorScaleStartTile };
    } else {
        $colorScaleStartTile = 0.025;
    }
    
    if( exists($opts->{ colorScaleEndTile }) ) {
        $colorScaleEndTile = $opts->{ colorScaleEndTile };
    } else {
        $colorScaleEndTile = 0.975;
    }
    
    if( exists($opts->{ posColorString }) ) {
        $posColorString = $opts->{ posColorString };
    } else {
        $posColorString="white,orange,red,darkRed";
    }    
    
    if( exists($opts->{ negColorString }) ) {
        $negColorString = $opts->{ negColorString };
    } else {
        $negColorString="white,cyan,blue,darkBlue";        
    }    
    
    if( exists($opts->{ missingColor }) ) {
        $missingColor = $opts->{ missingColor };
    } else {
        $missingColor="null";        
    }    
    
    if( exists($opts->{ highlightColor }) ) {
        $highlightColor = $opts->{ highlightColor };
    } else {
        $highlightColor="cyan";
    }
    
    if( exists($opts->{ elementBedFile }) ) {
        $elementBedFile = $opts->{ elementBedFile };
    } else {
        $elementBedFile="";
    }
    
    if( exists($opts->{ imageQuality }) ) {
        # 0 - 9 (0=best,9=worst)
        $imageQuality = $opts->{ imageQuality };
        $imageQuality = 9 if(($imageQuality < 0) or ($imageQuality > 9));
    } else {
        $imageQuality = 9
    }
    
    if( exists($opts->{ scaleFragmentSizes }) ) {
        $scaleFragmentSizes = 1;
    } else {
        $scaleFragmentSizes = 0;
    }
    
    if( exists($opts->{ drawTriangle }) ) {
        $drawTriangle = 1;
    } else {
        $drawTriangle = 0;
    }
    
    if( exists($opts->{ drawDiamond }) ) {
        $drawDiamond= 1;
    } else {
        $drawDiamond = 0;
    }
    
    if( exists($opts->{ transparentBGFlag }) ) {
        $transparentBGFlag = 1;
    } else {
        $transparentBGFlag = 0;
    }
    
    if( exists($opts->{ embed_meta }) ) {
        $embed_meta = 1;
    } else {
        $embed_meta = 0;
    }
    
    if( exists($opts->{ contigSpacing }) ) {
        $contigSpacing = $opts->{ contigSpacing };
        $contigSpacing = 1 if($contigSpacing < 1);
    } else {
        $contigSpacing = 1;
    }
    
    if( exists($opts->{ scaleMode }) ) {
        $scaleMode = $opts->{ scaleMode };
        $scaleMode = "combined" if(($scaleMode ne "combined") and ($scaleMode ne "seperate"));
    } else {
        $scaleMode="combined";
    }
    
    if( exists($opts->{ transparency }) ) {
        $transparency = $opts->{ transparency };
        $transparency = 0 if(($transparency < 0) || ($transparency > 255));
    } else {
        $transparency=0;
    }
    
    if(($drawTriangle) and ($drawDiamond)) {
        $drawTriangle=0;
        $drawDiamond=0;
        print STDERR "WARNING - cannot select both triangle and diamond options!\n";
        help();
    }
    
    $ret->{ inputMatrixArray }=$inputMatrixArray;
    $ret->{ verbose }=$verbose;
    $ret->{ output }=$output;
    $ret->{ outputSuffix }=$outputSuffix;
    $ret->{ imageSize }=$imageSize;
    $ret->{ pixelSize }=$pixelSize;
    $ret->{ y_pixelSize }=$y_pixelSize;
    $ret->{ x_pixelSize }=$x_pixelSize;
    $ret->{ maxImageDim }=$maxImageDim;
    $ret->{ drawPixelBorder }=$drawPixelBorder;
    $ret->{ omitContigBorder }=$omitContigBorder;
    $ret->{ drawLabel }=$drawLabel;
    $ret->{ drawScores }=$drawScores;
    $ret->{ logTransform }=$logTransform;
    $ret->{ colorScaleStart }=$colorScaleStart;
    $ret->{ colorScaleEnd }=$colorScaleEnd;
    $ret->{ colorScaleStartTile }=$colorScaleStartTile;
    $ret->{ colorScaleEndTile }=$colorScaleEndTile;
    $ret->{ posColorString }=$posColorString;
    $ret->{ negColorString }=$negColorString;
    $ret->{ missingColor }=$missingColor;
    $ret->{ highlightColor }=$highlightColor;
    $ret->{ elementBedFile }=$elementBedFile;
    $ret->{ imageQuality }=$imageQuality;
    $ret->{ scaleFragmentSizes }=$scaleFragmentSizes;
    $ret->{ drawTriangle }=$drawTriangle;
    $ret->{ drawDiamond }=$drawDiamond;
    $ret->{ transparentBGFlag }=$transparentBGFlag;
    $ret->{ embed_meta }=$embed_meta;
    $ret->{ contigSpacing }=$contigSpacing;
    $ret->{ scaleMode }=$scaleMode;
    $ret->{ transparency }=$transparency;
    
    return($ret,$inputMatrixArray,$verbose,$output,$outputSuffix,$imageSize,$pixelSize,$y_pixelSize,$x_pixelSize,$maxImageDim,$drawPixelBorder,$omitContigBorder,$drawLabel,$drawScores,$logTransform,$colorScaleStart,$colorScaleEnd,$colorScaleStartTile,$colorScaleEndTile,$posColorString,$negColorString,$missingColor,$highlightColor,$elementBedFile,$imageQuality,$scaleFragmentSizes,$drawTriangle,$drawDiamond,$transparentBGFlag,$embed_meta,$contigSpacing,$scaleMode,$transparency);
}

sub intro() {
    print STDERR "\n";
    
    print STDERR "Tool:\t\t".$tool."\n";
    print STDERR "Version:\t".$cworld::dekker::VERSION."\n";
    print STDERR "Summary:\tdraws heatmap PNG of matrix file\n";
    
    print STDERR "\n";
}

sub help() {
    intro();
    
    print STDERR "Usage: perl heatmap.pl [OPTIONS] -i <inputMatrix>\n";
    
    print STDERR "\n";
    
    print STDERR "Required:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-i@", "[]", "input matrix file [2 allowed, 1st=bottom-left, 2nd=top-right]");
    
    print STDERR "\n";
    
    print STDERR "Options:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-v", "[]", "FLAG, verbose mode");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-o", "[]", "prefix for output file(s)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--os", "[]", "optional output suffix for heatmap image");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--sfs", "[]", "FLAG, scaleFragmentSizes, scale pixels by fragment/bin size");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--dt", "[]", "FLAG, drawTriangle, draw heatmap as triangle (only upper triangle)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--dd", "[]", "FLAG, drawDiamond, draw heatmap as diamond");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--bg", "[]", "FLAG, use transparent background");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--em", "[]", "FLAG, embed meta data into PNG");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--dpb", "[]", "FLAG, draw pixel border");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ocb", "[]", "FLAG, omit contig/region border");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--dl", "[]", "FLAG, draw header labels");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ds", "[]", "FLAG, draw scores in pixels");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--name", "[]", "prefix for output file");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--is", "[800]", "ideal image size in pixel (larger of width/height");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ps", "[]", "x/y axis pixel size");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--yps", "[]", "y axis pixel size");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--xps", "[]", "x axis pixel size");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--maxdim", "[30000]", "maximum image dimension [hard-limited:32000]");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--lt", "[0]", "log transform input data into specified base");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--start", "[]", "absolute value for color start");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--end", "[]", "absolute vlaue for color end");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--startTile", "[0.025]", "fraction value for color start");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--endTile", "[0.975]", "fraction value for color end");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ebf", "[]", "element bed file - overlap row/col will be highlighted");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--iq", "[9]", "image quality, 0=best 9=worst");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--sm", "[]", "scale mode, combined or seperate for cis/trans [combined/seperate]");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--pc", "[white,orange,red,darkRed]", "positive color scale string");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--nc", "[white,cyan,blue,darkBlue]",  "negative color scale string");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--mc", "[null]", "missing data color");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--hc", "[cyan]", "highlight row/col data color");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--cs", "[]", "contig spacing, pixel size for spacing between contigs (contig border)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--t", "[]", "transparency setting for all colors [0-255] 0=none / 255=full");
    
    
    print STDERR "\n";
    
    print STDERR "Notes:";
    print STDERR "
    This script coverts a matrix into a PNG heatmap image.
    Matrix should have row/col headers in the standard my5C format.
    Matrix can be TXT or gzipped TXT.
    See website for matrix format details.\n";
    
    print STDERR "\n";
    
    print STDERR "Contact:
    Dekker Lab
    http://my5C.umassmed.edu
    my5C.help\@umassmed.edu\n";
    
    print STDERR "\n";
    
    exit;
}

sub getMaxScoreWidth($$) {
    my $matrixObject=shift;
    my $matrix=shift;
    
    my $inc2header=$matrixObject->{ inc2header };
    my $numYHeaders=$matrixObject->{ numYHeaders };
    my $numXHeaders=$matrixObject->{ numXHeaders };
    my $symmetrical = $matrixObject->{ symmetrical };
    my $missingValue = $matrixObject->{ missingValue };
    
    my $maxScoreWidth=0;
    for(my $y=0;$y<$numYHeaders;$y++) {
        my $yHeader=$inc2header->{ y }->{$y};
        for(my $x=0;$x<$numXHeaders;$x++) {
            my $xHeader=$inc2header->{ x }->{$x};
            
            # skip below diagonal
            next if(($symmetrical) and ($y > $x));
            
            my $score = $missingValue;
            $score=$matrix->{$y}->{$x} if(defined($matrix->{$y}->{$x}));
            
            next if($score eq "NA");
            next if( ($score eq $missingValue) or (($score ne "NA") and ($missingValue ne "NA") and ($score == $missingValue)) );
            
            $maxScoreWidth=length($score) if(length($score) > $maxScoreWidth);
            
        }
    }
    
    return($maxScoreWidth);
}

sub scaleByFragmentSize($$;$) {
    # required
    my $matrixObject=shift;
    my $pixelSize=shift;
    # optional
    my $headerSizeAccuracy=0.02;
    $headerSizeAccuracy = shift if @_;
    
    my $inc2header=$matrixObject->{ inc2header };
    my $header2inc=$matrixObject->{ header2inc };
    my $numYHeaders=$matrixObject->{ numYHeaders };
    my $numXHeaders=$matrixObject->{ numXHeaders };
    my $numTotalHeaders=$matrixObject->{ numTotalHeaders };
    my $verbose=$matrixObject->{ verbose };
    
    my %header2pixelSize=();
    
    # find min fragment size for X & Y axis
    my @headerSizes=();
    for(my $xy=0;$xy<$numTotalHeaders;$xy++) {
        my $header=$inc2header->{ xy }->{$xy};
        
        my $headerObject=getHeaderObject($header);
        my $headerSize=$headerObject->{ size };
    
        push(@headerSizes,$headerSize);
    }
    
    my $headerSizeStats=listStats(\@headerSizes,$headerSizeAccuracy);
    my $minHeaderSize=$headerSizeStats->{ trimmedMin };
    my $maxHeaderSize=$headerSizeStats->{ trimmedMax };
    
    print STDERR "\tminHeaderSize\t$minHeaderSize\n" if($verbose);
    print STDERR "\tmaxHeaderSize\t$maxHeaderSize\n" if($verbose);
    
    # get x axis header scaling factors
    my $imageWidth=0;
    for(my $x=0;$x<$numXHeaders;$x++) {
        my $header=$inc2header->{ x }->{$x};
        
        my $headerObject=getHeaderObject($header);
        my $headerSize=$headerObject->{ size };
    
        # handle bounds
        $headerSize=$minHeaderSize if($headerSize < $minHeaderSize);
        $headerSize=$maxHeaderSize if($headerSize > $maxHeaderSize);
        
        my $scaledPixelSize=ceil($headerSize/$minHeaderSize)*($pixelSize);
        
        #print STDERR "$header -> $headerSize / $minHeaderSize = $scaledPixelSize [$imageWidth]\n";
        $header2pixelSize{$header}=$scaledPixelSize;
        $imageWidth += $scaledPixelSize;
    }
    
    # get y axis header scaling factors
    my $imageHeight=0;
    for(my $y=0;$y<$numYHeaders;$y++) {
        my $header=$inc2header->{ y }->{$y};
        
        my $headerObject=getHeaderObject($header);
        my $headerSize=$headerObject->{ size };
        
        # handle bounds
        $headerSize=$minHeaderSize if($headerSize < $minHeaderSize);
        $headerSize=$maxHeaderSize if($headerSize > $maxHeaderSize);
    
        my $scaledPixelSize=ceil($headerSize/$minHeaderSize)*($pixelSize);
        
        #print STDERR "$header -> $headerSize / $minHeaderSize = $scaledPixelSize [$imageHeight]\n";
        $header2pixelSize{$header}=$scaledPixelSize;
        $imageHeight += $scaledPixelSize;
    }
        
    return(\%header2pixelSize,$imageWidth,$imageHeight);
    
}

my %options;
my $results = GetOptions( \%options,'inputMatrixArray|i=s@','verbose|v','output|o=s','outputSuffix|os=s','imageSize|is=i','pixelSize|ps=i','y_pixelSize|yps=i','x_pixelSize|xps=i','maxImageDim|maxdim=i','drawPixelBorder|dpb','omitContigBorder|ocb','drawLabel|dl','drawScores|ds','logTransform|lt=f','colorScaleStart|start=f','colorScaleEnd|end=f','colorScaleStartTile|startTile=f','colorScaleEndTile|endTile=f','posColorString|pc=s','negColorString|nc=s','missingColor|mc=s','highlightColor|hc=s','elementBedFile|ebf=s','imageQuality|iq=i','scaleFragmentSizes|sfs','drawTriangle|dt','drawDiamond|dd','transparentBGFlag|bg','embed_meta|em','contigSpacing|cs=f','scaleMode|sm=s','transparency|t=i') or croak help();
my ($ret,$inputMatrixArray,$verbose,$output,$outputSuffix,$imageSize,$pixelSize,$y_pixelSize,$x_pixelSize,$maxImageDim,$drawPixelBorder,$omitContigBorder,$drawLabel,$drawScores,$logTransform,$colorScaleStart,$colorScaleEnd,$colorScaleStartTile,$colorScaleEndTile,$posColorString,$negColorString,$missingColor,$highlightColor,$elementBedFile,$imageQuality,$scaleFragmentSizes,$drawTriangle,$drawDiamond,$transparentBGFlag,$embed_meta,$contigSpacing,$scaleMode,$transparency)=check_options( \%options );

intro() if($verbose);

#get the absolute path of the script being used
my $cwd = getcwd();
my $fullScriptPath=abs_path($0);
my @fullScriptPathArr=split(/\//,$fullScriptPath);
@fullScriptPathArr=@fullScriptPathArr[0..@fullScriptPathArr-3];
my $scriptPath=join("/",@fullScriptPathArr);
my $commentLine=getScriptOpts($ret,$tool);

my ($inputMatrix);
my ($matrix);
my ($matrixObject);

# for double heatmap
if(@{$inputMatrixArray} == 2) {
    
    my ($inputMatrix_1,$inputMatrix_2)=@{$inputMatrixArray};
    
    $output="default" if($output eq "");
    $output .= "__".$outputSuffix if($outputSuffix ne "");
    
    # validate all files are the same
    print STDERR "validating identical matrices...\n" if($verbose);
    for(my $i1=0;$i1<@{$inputMatrixArray};$i1++) {
        my $inputMatrix_1 = $inputMatrixArray->[$i1];
        
        croak "inputMatrix [$inputMatrix_1] does not exist" if(!(-e $inputMatrix_1));
        
        for(my $i2=0;$i2<@{$inputMatrixArray};$i2++) {
            my $inputMatrix_2 = $inputMatrixArray->[$i2];
            
            next if($i1 >= $i2);
            croak "inputMatrix [$inputMatrix_2] does not exist" if(!(-e $inputMatrix_2));
            
            validateIdenticalMatrixStructure($inputMatrix_1,$inputMatrix_2);
        }
    }
    print STDERR "\tdone\n" if($verbose);
    
    croak "inputMatrix_1 must be symmetrical! [$inputMatrix_1]" if(!isSymmetrical($inputMatrix_1));
    croak "inputMatrix_2 must be symmetrical! [$inputMatrix_2]" if(!isSymmetrical($inputMatrix_2));

    print STDERR "\n" if($verbose);
    
    # get matrix 1 information
    my $matrixObject_1=getMatrixObject($inputMatrix_1,$output,$verbose);
    my $inc2header_1=$matrixObject_1->{ inc2header };
    my $header2inc_1=$matrixObject_1->{ header2inc };
    my $numYHeaders_1=$matrixObject_1->{ numYHeaders };
    my $numXHeaders_1=$matrixObject_1->{ numXHeaders };
    
    my ($matrix_1)=getData($inputMatrix_1,$matrixObject_1,$verbose);
    
    if($logTransform > 0) {
        print STDERR "\nlog transformming matrix (lt=$logTransform)...\n" if($verbose);
        $matrix_1=logTransformMatrix($matrix_1,$matrixObject_1,$logTransform);
        print STDERR "\tdone\n" if($verbose);
    }

    print STDERR "\n" if($verbose);

    # get matrix 2 information
    my $matrixObject_2=getMatrixObject($inputMatrix_2,$output,$verbose);
    my $inc2header_2=$matrixObject_2->{ inc2header };
    my $header2inc_2=$matrixObject_2->{ header2inc };
    my $numYHeaders_2=$matrixObject_2->{ numYHeaders };
    my $numXHeaders_2=$matrixObject_2->{ numXHeaders };

    my ($matrix_2)=getData($inputMatrix_2,$matrixObject_2,$verbose);

    if($logTransform > 0) {
        print STDERR "\nlog transformming matrix (lt=$logTransform)...\n" if($verbose);
        $matrix_2=logTransformMatrix($matrix_2,$matrixObject_2,$logTransform);
        print STDERR "\tdone\n" if($verbose);
    }

    print STDERR "\n" if($verbose);

    my $inputMatrixName_1=getFileName($inputMatrix_1);
    my $inputMatrixName_2=getFileName($inputMatrix_2);
    my $inputMatrixName = $inputMatrixName_1."___".$inputMatrixName_2;
    $inputMatrixName=$output if($output ne "");

    # stitch matrices
    print STDERR "stiching matrices ...\n" if($verbose);
    ($matrix,my $stitch_missingValue)=stitchMatrices($matrixObject_1,$matrixObject_2,$matrix_1,$matrix_2);
    print STDERR "\t$stitch_missingValue\n" if($verbose);
    print STDERR "\tdone\n" if($verbose);
    

    $matrix_1=undef;
    $matrix_2=undef;
    
    print STDERR "\n" if($verbose);

    my $stitchMatrixFile=$output.".double.matrix.gz";
    print STDERR "writing stitched matrix ($stitchMatrixFile)...\n" if($verbose);
    writeMatrix($matrix,$inc2header_1,$stitchMatrixFile,$stitch_missingValue,$commentLine);
    print STDERR "\tdone\n" if($verbose);
    
    print STDERR "\n" if($verbose);
    
    $inputMatrix=$stitchMatrixFile;
    
    $matrixObject=getMatrixObject($inputMatrix,$output,$verbose);
    
} elsif(@{$inputMatrixArray}) {
    
    $inputMatrix=$inputMatrixArray->[0];
    
    $matrixObject=getMatrixObject($inputMatrix,$output,$verbose);
    
    ($matrix)=getData($inputMatrix,$matrixObject,$verbose);

} else {
    
   croak "must supply 1 or 2 matrices";

}

print STDERR "\n" if($verbose);

my $inc2header=$matrixObject->{ inc2header };
my $header2inc=$matrixObject->{ header2inc };
my $numYHeaders=$matrixObject->{ numYHeaders };
my $numXHeaders=$matrixObject->{ numXHeaders };
my $numTotalHeaders=$matrixObject->{ numTotalHeaders };
my $yHeaderLength=$matrixObject->{ yHeaderLength };
my $xHeaderLength=$matrixObject->{ xHeaderLength };
my $missingValue=$matrixObject->{ missingValue };
my $symmetrical=$matrixObject->{ symmetrical };
my $nInteractions=($numYHeaders*$numXHeaders);
my $numContigs=$matrixObject->{ numContigs };
my $numXContigs=$matrixObject->{ numXContigs };
my $numYContigs=$matrixObject->{ numYContigs };
my $inputMatrixName=$matrixObject->{ inputMatrixName };
$output=$matrixObject->{ output };

$output .= "__".$outputSuffix if($outputSuffix ne "");

if($logTransform > 0) {
    print STDERR "log transformming matrix (lt=$logTransform)...\n" if($verbose);
    $matrix=logTransformMatrix($matrix,$matrixObject,$logTransform);
    print STDERR "\tdone\n" if($verbose);
    
    print STDERR "\n" if($verbose);

}

# calculate the color scale
my ($cisColorScaleStart,$cisColorScaleEnd,$transColorScaleStart,$transColorScaleEnd);
$cisColorScaleStart=$transColorScaleStart=$colorScaleStart;
$cisColorScaleEnd=$transColorScaleEnd=$colorScaleEnd;

if(($colorScaleStart eq "NA") or ($colorScaleEnd eq "NA")) {
    print STDERR "Calculating Color Scale [$colorScaleStartTile:$colorScaleEndTile:$scaleMode] ...\n" if($verbose);
    if($scaleMode eq "seperate") {
        ($cisColorScaleStart,$cisColorScaleEnd)=autoScale($matrixObject,$matrix,$colorScaleStartTile,$colorScaleEndTile,"cis");
        print STDERR "\tCIS\t".$cisColorScaleStart."-".$cisColorScaleEnd."\n" if($verbose);
        ($transColorScaleStart,$transColorScaleEnd)=autoScale($matrixObject,$matrix,$colorScaleStartTile,$colorScaleEndTile,"trans");
        print STDERR "\tTRANS\t".$transColorScaleStart."-".$transColorScaleEnd."\n" if($verbose);
    } else {
        my ($tmpColorScaleStart,$tmpColorScaleEnd)=autoScale($matrixObject,$matrix,$colorScaleStartTile,$colorScaleEndTile,"all");
        $colorScaleStart=$tmpColorScaleStart if($colorScaleStart eq "NA");
        $colorScaleEnd=$tmpColorScaleEnd if($colorScaleEnd eq "NA");
        print STDERR "\tALL\t".$colorScaleStart."-".$colorScaleEnd."\n" if($verbose);
    }
    print STDERR "\n" if($verbose);
}

my $textOffset=5;
my $labelFont=GD::Font->Large;
my ($labelFontWidth,$labelFontHeight) = ($labelFont->width, $labelFont->height);
my $scoreFont=GD::Font->Tiny;
my ($scoreFontWidth,$scoreFontHeight) = ($scoreFont->width, $scoreFont->height);

# setup image parameters
$imageSize=$maxImageDim if($imageSize > $maxImageDim);
my $colorString=$posColorString."___".$negColorString;
my $labelPixelSize=max($labelFontWidth,$labelFontHeight);

#calculate auto pixel size
print STDERR "Calculating Auto Pixel Size ($imageSize)...\n" if($verbose);
$pixelSize=autoSize($imageSize,$numYHeaders,$numXHeaders) if($pixelSize eq "NA");
# disable draw label if image_size > max_dim
if( ($drawLabel) and ((($numYHeaders * $labelPixelSize) > $maxImageDim) or  (($numXHeaders * $labelPixelSize) > $maxImageDim) ) ) {
    print STDERR "\tCannot enable row/col labels [".$drawLabel."] - image too large (".$pixelSize."->".$labelPixelSize.") [".($numYHeaders*$labelPixelSize).">".$maxImageDim." or ".($numXHeaders*$labelPixelSize).">".$maxImageDim."]\n";
    $drawLabel=0;
}
$pixelSize=$labelPixelSize if(($drawLabel) and (($numYHeaders * $labelPixelSize) < $maxImageDim) and (($numXHeaders * $labelPixelSize) < $maxImageDim) );

my $maxScoreWidth=0;
($maxScoreWidth)=getMaxScoreWidth($matrixObject,$matrix) if($drawScores);
$pixelSize = (($scoreFontWidth*$maxScoreWidth)+($textOffset*2)) if($drawScores);
$y_pixelSize=$pixelSize if($y_pixelSize eq "NA");
$x_pixelSize=$pixelSize if($x_pixelSize eq "NA");

my $imageHeight=($numYHeaders*$y_pixelSize);
my $imageWidth=($numXHeaders*$x_pixelSize);
croak "heatmap too large! [$imageHeight x $imageWidth > $maxImageDim]" if(($imageHeight > $maxImageDim) or ($imageWidth > $maxImageDim));
print STDERR "\tpixelSize [p|y|x] [".$pixelSize."|".$y_pixelSize."|".$x_pixelSize."] ($imageHeight x $imageWidth)\n" if($verbose);

print STDERR "\n" if($verbose);

# process optional fragment scale flag (scale pixels by fragment length)
my ($header2pixelSize);
my $headerSizeAccuracy=0.02; # init with 2% trimming
if($scaleFragmentSizes) {
    print STDERR "Scaling pixel size by fragment size ...\n" if($verbose);
    
    # reset pixel size
    $pixelSize=$y_pixelSize=$x_pixelSize=1;
    
    ($header2pixelSize,$imageWidth,$imageHeight)=scaleByFragmentSize($matrixObject,$pixelSize,$headerSizeAccuracy);
    
    while((($imageHeight > $maxImageDim) or ($imageWidth > $maxImageDim)) and ($headerSizeAccuracy < 1)) {
        $headerSizeAccuracy += 0.01;
        print "\timage too large - descreasing accuracy ... \n" if($verbose);
        ($header2pixelSize,$imageWidth,$imageHeight)=scaleByFragmentSize($matrixObject,$pixelSize,$headerSizeAccuracy);
    }
    
    croak "heatmap too large! [beyond maxImageDimxmaxImageDim PNG limit]" if(($imageHeight > $maxImageDim) or ($imageWidth > $maxImageDim));
    
    print STDERR "\n" if($verbose);
    $output .= ".fragmentScaled" if($scaleFragmentSizes);
}

$contigSpacing += 1 if( ($contigSpacing != 0) and (($contigSpacing % 2) == 0) and ($omitContigBorder == 0) and ($numContigs > 1) );

# fine tune image dimensions (pixel border | pixel label)
print STDERR "Initializing Image...\n" if($verbose);
print STDERR "\tdrawPixelBorder [$drawPixelBorder]\n" if($verbose);
print STDERR "\tomitContigBorder [$omitContigBorder][$contigSpacing]\n" if($verbose);
print STDERR "\t[p|y|x] [".$pixelSize."|".$y_pixelSize."|".$x_pixelSize."]\n" if($verbose);
print STDERR "\timageDimensions=($imageHeight x $imageWidth)\n" if($verbose);
# increase size for pixel borders
$imageHeight += ($numYHeaders*2) if($drawPixelBorder);
$imageWidth += ($numXHeaders*2) if($drawPixelBorder);
print STDERR "\tpixelBorder\textendedImageDimensions=($imageHeight x $imageWidth)\n" if($verbose);
# increase size for contig borders
$imageHeight += ($numYContigs)*($contigSpacing)+1 if(($omitContigBorder == 0) and ($numContigs > 1));
$imageWidth += ($numXContigs)*($contigSpacing)+1 if(($omitContigBorder == 0) and ($numContigs > 1));
print STDERR "\tcontigBorder\textendedImageDimensions=($imageHeight x $imageWidth)\n" if($verbose);
# log matrix dimensions
my $matrixHeight=$imageHeight;
my $matrixWidth=$imageWidth;
# increase size for labels
$textOffset += 1 if($drawPixelBorder);
$textOffset += ceil($contigSpacing/2) if($drawPixelBorder);
my $maxTextHeight=($xHeaderLength*$labelFontWidth)+($textOffset*2);
my $maxTextWidth=($yHeaderLength*$labelFontWidth)+($textOffset*2);
$imageHeight += $maxTextHeight if($drawLabel);
$imageWidth += $maxTextWidth if($drawLabel);
print STDERR "\trowcol labels\textendedImageDimensions=($imageHeight x $imageWidth)\n" if($verbose);
print STDERR "\n" if($verbose);

my $highlightHeaders={};
if(-e($elementBedFile)) {
    print STDERR "Intersecting headers with highlight bed file...\n" if($verbose);
    $highlightHeaders=intersectHeaders($matrixObject,$elementBedFile);
    print STDERR "\t" . keys(%{$highlightHeaders}) . " highlight headers\n" if($verbose);
    print STDERR "\n" if($verbose);
}

# create the perl gd object
print STDERR "Calculating Color Information...\n" if($verbose);
my ($img,$colorPalette,$nColorShades,$availableColors)=initHeatmap($imageHeight,$imageWidth,$colorString,$missingColor,$highlightColor,$transparency,$verbose);

# calculate color distances
my ($colorDistance,$colorBucketSize,$cisColorDistance,$cisColorBucketSize,$transColorDistance,$transColorBucketSize);
$colorDistance=$colorBucketSize=$cisColorDistance=$cisColorBucketSize=$transColorDistance=$transColorBucketSize="NA";
$colorDistance = ($colorScaleEnd-$colorScaleStart) if(($colorScaleEnd ne "NA") and ($colorScaleStart ne "NA"));
$colorBucketSize = ($colorDistance/($nColorShades-1)) if($colorDistance ne "NA");
$cisColorDistance = ($cisColorScaleEnd-$cisColorScaleStart) if(($cisColorScaleEnd ne "NA") and ($cisColorScaleStart ne "NA"));
$cisColorBucketSize = ($cisColorDistance/($nColorShades-1)) if($cisColorDistance ne "NA");
$transColorDistance = ($transColorScaleEnd-$transColorScaleStart) if(($transColorScaleEnd ne "NA") and ($transColorScaleStart ne "NA"));
$transColorBucketSize = ($transColorDistance/($nColorShades-1)) if($transColorDistance ne "NA");
print STDERR "\t$colorDistance [$colorBucketSize]\t$cisColorDistance [$cisColorBucketSize]\t$transColorDistance [$transColorBucketSize]\n" if($verbose);
print STDERR "\n" if($verbose);

# init GD brush positions
my $brushX=0;
$brushX=1 if($drawPixelBorder);
my $brushY=0;
$brushY=$maxTextHeight if($drawLabel);
$brushY += 1 if($drawPixelBorder);

# log matrix bounds
my $matrixStart_y=$brushY;
my $matrixStart_x=$brushX;

# x/y offset for pixel border
my $yOffset=0;
$yOffset=1 if($drawPixelBorder);
my $xOffset=0;
$xOffset=1 if($drawPixelBorder);

# init last regioins
my $lastYRegion="NA";
my $lastXRegion="NA";

# init progress meter
my $progressBucketSize=ceil($numYHeaders / 1000);
my $pcComplete=0;

my $headerObjects={};

my $drawMode="rectangle";
$drawMode="pixel" if(($x_pixelSize == 0) and ($y_pixelSize == 0));

# draw heatmap
print STDERR "drawing heatmap image [$drawMode] ...\n" if($verbose);
my ($y,$x);
for($y=0;$y<$numYHeaders;$y++) {
    $brushX=0;
    $brushX=1 if($drawPixelBorder);
    
    my $yHeader=$inc2header->{ y }->{$y};
    my $yHeaderObject={};
    if(!exists($headerObjects->{$yHeader})) {
        $yHeaderObject=getHeaderObject($yHeader);
        $headerObjects->{$yHeader}=$yHeaderObject;
    } else {
        $yHeaderObject=$headerObjects->{$yHeader};
    }
    my $yHeaderRegion=$yHeaderObject->{ region };
    
    my $tmp_y_pixelSize = $y_pixelSize;
    $tmp_y_pixelSize = $header2pixelSize->{$yHeader} if(($scaleFragmentSizes) and (exists($header2pixelSize->{$yHeader})));
    
    # draw region boundaries
    if(($omitContigBorder == 0) and ($lastYRegion ne $yHeaderRegion) and ($numContigs > 1)) {
        my $brush_regionY = $brushY;
        $brush_regionY = $brushY + floor($contigSpacing/2) - $yOffset if($lastYRegion ne "NA");
        $brush_regionY -= $yOffset if($lastYRegion eq "NA");
        $brushY+=$contigSpacing if($lastYRegion ne "NA");
        $brushY+=ceil($contigSpacing/2) if($lastYRegion eq "NA");
        $img->line(($matrixStart_x-$xOffset),$brush_regionY,(($matrixStart_x-$xOffset)+$matrixWidth)-1,$brush_regionY,$colorPalette->{ B });
    }
        
    $lastXRegion="NA";

    for($x=0;$x<$numXHeaders;$x++) {
    
        my $xHeader=$inc2header->{ x }->{$x}; 
        my $xHeaderObject={};
        if(!exists($headerObjects->{$xHeader})) {
            $xHeaderObject=getHeaderObject($xHeader);
            $headerObjects->{$xHeader}=$xHeaderObject;
        } else {
            $xHeaderObject=$headerObjects->{$xHeader};
        }
        my $xHeaderRegion=$xHeaderObject->{ region };
        
        my $tmp_x_pixelSize = $x_pixelSize;
        $tmp_x_pixelSize = $header2pixelSize->{$xHeader} if(($scaleFragmentSizes) and (exists($header2pixelSize->{$xHeader})));
        
        my $score=$matrixObject->{ missingValue };
        $score=$matrix->{$y}->{$x} if(defined($matrix->{$y}->{$x}));
        
        my $colorIndex = -1;
        $colorIndex=getColorIndex($score,$colorScaleStart,$colorScaleEnd,$nColorShades,$colorBucketSize) if(($score ne "NA") and ($scaleMode eq "combined"));
        if($scaleMode eq "seperate") {
            my $interactionDistance=getInteractionDistance($matrixObject,$yHeaderObject,$xHeaderObject,1);
            $colorIndex = getColorIndex($score,$cisColorScaleStart,$cisColorScaleEnd,$nColorShades,$cisColorBucketSize) if(($interactionDistance != -1) and ($score ne "NA"));
            $colorIndex = getColorIndex($score,$transColorScaleStart,$transColorScaleEnd,$nColorShades,$transColorBucketSize) if(($interactionDistance == -1) and ($score ne "NA"));
        }
            
        my $color=$colorPalette->{ N };
        
        if($score eq "NA") {
            $color=$colorPalette->{ NA } if(exists($colorPalette->{ NA }));
        } elsif($score < 0) {
            $color=$colorPalette->{ nc }[$colorIndex] if(($colorIndex != -1) and (exists($colorPalette->{ nc }[$colorIndex])));
        } elsif($score >= 0) {
            $color=$colorPalette->{ pc }[$colorIndex] if(($colorIndex != -1) and (exists($colorPalette->{ pc }[$colorIndex])));
        }
        
        # draw region boundaries
        if(($omitContigBorder == 0) and ($lastXRegion ne $xHeaderRegion) and ($numContigs > 1)) {
            my $brush_regionX = $brushX;
            $brush_regionX = $brushX + floor($contigSpacing/2) - $xOffset if($lastXRegion ne "NA");
            $brush_regionX -= $xOffset if($lastXRegion eq "NA");
            $brushX+=$contigSpacing if($lastXRegion ne "NA");
            $brushX+=ceil($contigSpacing/2) if($lastXRegion eq "NA");
            $img->line($brush_regionX,($matrixStart_y-$yOffset),$brush_regionX,(($matrixStart_y-$yOffset)+$matrixHeight)-1,$colorPalette->{ B }) if($y == ($numYHeaders-1));
        }
        
        #draw pixels
        $img->setPixel($brushX,$brushY,$color) if($drawMode eq "pixel");
        $img->filledRectangle(($brushX),($brushY),($brushX+$tmp_x_pixelSize),($brushY+$tmp_y_pixelSize),$color) if($drawMode eq "rectangle");
        
        # optional row/col highlight
        $img->filledRectangle(($brushX),($brushY),($brushX+$tmp_x_pixelSize),($brushY+$tmp_y_pixelSize),$colorPalette->{ HIGHLIGHT }) if( ((exists($highlightHeaders->{$yHeader})) or (exists($highlightHeaders->{$xHeader})) ) );
        
        #draw pixel border
        $img->rectangle($brushX-1,$brushY-1,($brushX+$tmp_x_pixelSize+1),($brushY+$tmp_y_pixelSize+1),$colorPalette->{ G }) if($drawPixelBorder);
        
        $img->string($scoreFont,($brushX+ceil($tmp_x_pixelSize-round(length($score)*$scoreFontWidth))/2),round($brushY+floor($tmp_y_pixelSize/2)),$score,$colorPalette->{ B }) if($drawScores);
        
        $brushX+=2 if($drawPixelBorder);
        
        $img->stringUp($labelFont,($brushX+($tmp_x_pixelSize/2))-($labelFontHeight/2),$brushY-$textOffset,$xHeader,$colorPalette->{ B }) if(($y == 0) and ($drawLabel));
        
        $lastXRegion=$xHeaderRegion;
        
        $brushX += $tmp_x_pixelSize;
    
    }
    
    $img->string($labelFont,$brushX+$textOffset,($brushY+($tmp_y_pixelSize/2))-($labelFontHeight/2),$yHeader,$colorPalette->{ B }) if($drawLabel);
    
    $brushY += $tmp_y_pixelSize;
    $brushY+=2 if($drawPixelBorder);
    
    $lastYRegion=$yHeaderRegion;
    
    $pcComplete = 100 if($y == ($numYHeaders-1));
    print STDERR "\e[A" if(($y != 0) and ($verbose));
    printf STDERR "\t%.2f%% complete ($y / ".($numYHeaders-1).")...\n", $pcComplete if($verbose);
    $pcComplete = round((($y/$numYHeaders)*100),2);

}

# draw region boundaries
if(($omitContigBorder == 0) and ($numContigs > 1)) {
    my $brush_regionY = 0;
    $brush_regionY=$brushY + floor($contigSpacing/2) - $yOffset;
    $img->line(($matrixStart_x-$xOffset),$brush_regionY,(($matrixStart_x-$xOffset)+$matrixWidth)-1,$brush_regionY,$colorPalette->{ B });
}
    
print STDERR "\tcomplete\n" if($verbose);

print STDERR "\n" if($verbose);

if( (($drawTriangle) or ($drawDiamond)) and (!$symmetrical)) {
    
    print STDERR "error - triangle/diamond requires a symmetrical matrix!\n";
    print STDERR "\n" if($verbose);    
    
} elsif( (($drawTriangle) or ($drawDiamond)) and ($symmetrical)) {

    $output .= ".triangle" if($drawTriangle);
    $output .= ".diamond" if($drawDiamond);

    # old trig method
    # calculate new rotated image dimensions
    #my $triangleArea=ceil(($imageHeight*$imageWidth)/2)
    #my $newImageWidth = ceil(sqrt( ($imageWidth**2) + ($imageHeight**2) ));
    #y $newImageHeight = ceil(2*($triangleArea/$newImageWidth));
    ##opp/hyp
    #my $triangleAngle=(asin($imageWidth/$newImageWidth)*(180/PI));
    #$triangleAngle=45; # hard coded for symmetrical maps

    #rotate
    my $diagonalSize = sqrt( ($pixelSize**2) + ($pixelSize**2) );
    my $newImageWidth = ceil( $imageWidth * ($diagonalSize/$pixelSize) );
    my $newImageHeight = ceil( $imageHeight * ($diagonalSize/$pixelSize) );
    $newImageHeight = ceil($newImageHeight/2) if($drawTriangle);
    
    $newImageHeight += 1 if(($newImageHeight % 2));
    $newImageWidth += 1 if(($newImageWidth % 2));
    my $newImg = new GD::Image($newImageWidth,$newImageHeight);
    
    $newImg->fill(0,0,$colorPalette->{ BG });
    
    # draw rotated diamond
    $newImg->copyRotated($img,($newImageWidth/2),($newImageHeight/2),0,0,$imageWidth,$imageHeight,45) if($drawDiamond);
    
    # draw rotated upper triangle
    $newImg->copyRotated($img,($newImageWidth/2),($newImageHeight),0,0,$imageWidth,$imageHeight,45) if($drawTriangle);
    
    $img=$newImg;
}

# hide background
if($transparentBGFlag) {
    print STDERR "Converting background color to transparent...\n" if($verbose);
    $img->transparent($colorPalette->{ BG });
    print STDERR "\tcomplete\n\n" if($verbose);
}

my $pngFile=$output.".png";
print STDERR "Converting object to png ($inputMatrixName)...\n" if($verbose);

open(OUT,">",$pngFile) or croak "Could not open file [$pngFile] - $!";
chmod(0644, $pngFile);
binmode OUT;
print OUT $img->png($imageQuality);
close(OUT);

print STDERR "\tcomplete\n\n" if($verbose);

if($embed_meta) {
    my $computeResource=getComputeResource();
    my $matrixSum=getMatrixSum($matrixObject,$matrix);
    
    system("convert '".$pngFile."' -set cw_inputMatrix '".$inputMatrix."' '".$pngFile."'");
    system("convert '".$pngFile."' -set cw_numYHeaders '".$numYHeaders."' '".$pngFile."'");
    system("convert '".$pngFile."' -set cw_numXHeader '".$numXHeaders."' '".$pngFile."'");
    system("convert '".$pngFile."' -set cw_pixelSize '".$pixelSize."' '".$pngFile."'");
    system("convert '".$pngFile."' -set cw_colorScaleStart '".$colorScaleStart."' '".$pngFile."'");
    system("convert '".$pngFile."' -set cw_colorScaleEnd '".$colorScaleEnd."' '".$pngFile."'");
    system("convert '".$pngFile."' -set cw_computeResource '".$computeResource."' '".$pngFile."'");
    system("convert '".$pngFile."' -set cw_missingValue '".$missingValue."' '".$pngFile."'");
    system("convert '".$pngFile."' -set cw_posColorString '".$posColorString."' '".$pngFile."'");
    system("convert '".$pngFile."' -set cw_negColorString '".$negColorString."' '".$pngFile."'");
    system("convert '".$pngFile."' -set cw_elementBedFile '".$elementBedFile."' '".$pngFile."'");
    system("convert '".$pngFile."' -set cw_highlightColor '".$highlightColor."' '".$pngFile."'");
    system("convert '".$pngFile."' -set cw_matrixSum '".$matrixSum."' '".$pngFile."'");
    
}