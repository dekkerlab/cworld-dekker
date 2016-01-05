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

use cworld::dekker;

my $tool=(split(/\//,abs_path($0)))[-1];

sub check_options {
    my $opts = shift;

    my ($inputMatrix,$verbose,$output,$minDistance,$maxDistance,$imageWidth,$headerSizingFile,$colorScaleStartTile,$colorScaleEndTile,$colorScaleStart,$colorScaleEnd,$posColorString,$negColorString,$missingColorString,$drawBGInteractions);
    
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
        $maxDistance = undef;
    }
    
    if( exists($opts->{ imageWidth }) ) {
        $imageWidth = $opts->{ imageWidth };
    } else {
        $imageWidth=1800;
    }
    
    if( exists($opts->{ headerSizingFile }) ) {
        $headerSizingFile = $opts->{ headerSizingFile };
    } else {
        $headerSizingFile = "";
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
    
    if( exists($opts->{ missingColorString }) ) {
        $missingColorString = $opts->{ missingColorString };
    } else {
        $missingColorString="N";
    }    
    
    if( exists($opts->{ drawBGInteractions }) ) {
        $drawBGInteractions = 1;
    } else {
        $drawBGInteractions = 0;
    }
    
    $ret->{ inputMatrix }=$inputMatrix;
    $ret->{ verbose }=$verbose;
    $ret->{ output }=$output;
    $ret->{ minDistance }=$minDistance;
    $ret->{ maxDistance }=$maxDistance;
    $ret->{ imageWidth }=$imageWidth;
    $ret->{ headerSizingFile }=$headerSizingFile;
    $ret->{ colorScaleStartTile }=$colorScaleStartTile;
    $ret->{ colorScaleEndTile }=$colorScaleEndTile;
    $ret->{ colorScaleStart }=$colorScaleStart;
    $ret->{ colorScaleEnd }=$colorScaleEnd;
    $ret->{ posColorString }=$posColorString;
    $ret->{ negColorString }=$negColorString;
    $ret->{ missingColorString }=$missingColorString;
    $ret->{ drawBGInteractions }=$drawBGInteractions;
    
    return($ret,$inputMatrix,$verbose,$output,$minDistance,$maxDistance,$imageWidth,$headerSizingFile,$colorScaleStartTile,$colorScaleEndTile,$colorScaleStart,$colorScaleEnd,$posColorString,$negColorString,$missingColorString,$drawBGInteractions);
}


sub intro() {
    print STDERR "\n";
    
    print STDERR "Tool:\t\t".$tool."\n";
    print STDERR "Version:\t".$cworld::dekker::VERSION."\n";
    print STDERR "Summary:\tdraws 'web-plot' of matrix file\n";
    
    print STDERR "\n";
}

sub help() {
    intro();
    
    print STDERR "Usage: perl matrix2webplot.pl [OPTIONS] -i <inputMatrix>\n";
    
    print STDERR "\n";
    
    print STDERR "Required:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-i@", "[]", "input matrix file [MULTIPLE]");
    
    print STDERR "\n";
    
    print STDERR "Options:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-v", "[]", "FLAG, verbose mode");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-o", "[]", "prefix for output file(s)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--minDist", "[]", "minimum allowed interaction distance, exclude < N distance (in BP)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--maxDist", "[]", "maximum allowed interaction distance, exclude > N distance (in BP)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--bg", "[]", "FLAG, use transparent background");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--maxdim", "[16000]", "maximum image dimension [hard-limited:32000]");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--lt", "[0]", "log transform input data into specified base");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--start", "[]", "absolute value for color start");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--end", "[]", "absolute vlaue for color end");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--startTile", "[0.025]", "fraction value for color start");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--endTile", "[0.975]", "fraction value for color end");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--pc", "[white,orange,red,darkRed]", "positive color scale string");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--nc", "[white,cyan,blue,darkBlue]",  "negative color scale string");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--mc", "[null]", "missing data color");
    
    
    print STDERR "\n";
    
    print STDERR "Notes:";
    print STDERR "
    This script coverts a matrix into a PNG matrix2webplot image.
    Matrix should have row/col headers in the standard my5C format.
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

sub initImage($$$$$) {
    my $colorScaleStart=shift;
    my $colorScaleEnd=shift;
    my $imageHeight=shift;
    my $imageWidth=shift;
    my $colorString=shift;
    
    GD::Image->trueColor(1);
    my $img = new GD::Image($imageWidth,$imageHeight);
    
    my %availableColors=();

    my $transparency=126;
    
    # utility colorPalette
    @{$availableColors{ null }} = (150,150,150,$transparency);
    @{$availableColors{ grey }} = (100,100,100,$transparency);
    @{$availableColors{ background }} = (225,225,225,120);
    #basic
    @{$availableColors{ white }} = (255,255,255,$transparency);
    @{$availableColors{ black }} = (0,0,0,$transparency);
    #primary
    @{$availableColors{ red }} = (255,0,0,$transparency);
    @{$availableColors{ green }} = (0,255,0,$transparency);
    @{$availableColors{ blue }} = (0,0,255,$transparency);
    #secondary
    @{$availableColors{ yellow }} = (255,255,0,$transparency);
    @{$availableColors{ violet }} = (255,0,255,$transparency);
    @{$availableColors{ cyan }} = (0,255,255,$transparency);
    #special 
    @{$availableColors{ orange }} = (255,127,0,$transparency);
    @{$availableColors{ blueGreen }} = (0,127,255,$transparency);
    #special 
    @{$availableColors{ yellowGreen }} = (0,255,127,$transparency);
    @{$availableColors{ redViolet }} = (255,0,127,$transparency);
    #special 
    @{$availableColors{ blueViolet }} = (127,0,255,$transparency);
    @{$availableColors{ yellowOrange }} = (127,255,0,$transparency);
    #special 
    @{$availableColors{ indigo }} = (75,0,130,$transparency);
    #special 
    @{$availableColors{ darkRed }} = (128,0,0,$transparency);
    @{$availableColors{ darkGreen }} = (0,128,0,$transparency);
    @{$availableColors{ darkBlue }} = (0,0,128,$transparency);
    #monoChromatic
    @{$availableColors{ monoRedStart }} = (255,203,203,$transparency);
    @{$availableColors{ monoRedEnd }} = (5,0,0,$transparency);
    @{$availableColors{ momoGreenStart }} = (233,255,203,$transparency);
    @{$availableColors{ monoGreenEnd }} = (0,5,0,$transparency);
    @{$availableColors{ monoBlueStart }} = (203,203,255,$transparency);
    @{$availableColors{ monoBlueEnd }} = (0,0,5,$transparency);
    
    my $colorPalette={};

    $colorPalette->{ W } = $img->colorAllocateAlpha(@{$availableColors{ white }});
    $colorPalette->{ N } = $img->colorAllocateAlpha(@{$availableColors{ null }});
    $colorPalette->{ G } = $img->colorAllocateAlpha(@{$availableColors{ grey }});
    $colorPalette->{ B } = $img->colorAllocateAlpha(@{$availableColors{ black }});
    $colorPalette->{ BG } = $img->colorAllocateAlpha(@{$availableColors{ background }});
    
    $img->fill(0,0,$colorPalette->{ W });
    
    my $nColorShades = 255;
    $nColorShades -= keys(%{$colorPalette});
    #print STDERR "\tnColorShades\t$nColorShades\n";
    
    my $nPositiveColorShades=floor($nColorShades/2);
    my $nNegativeColorShades=floor($nColorShades/2);
    #print STDERR "\tinitial color split\tnPosShades=$nPositiveColorShades\tnNegShades=$nNegativeColorShades\n";
    
    my ($positiveColorString,$negativeColorString)=split(/___/,$colorString);
    
    # positive color information
    $positiveColorString=$positiveColorString.",".$positiveColorString if($positiveColorString !~ /,/);
    my @positiveColorStringArray=split(/,/,$positiveColorString);
    my $nPositiveColors=@positiveColorStringArray;
    my $nPositiveShadesPerColor=floor($nPositiveColorShades/($nPositiveColors-1));    
    $nPositiveColorShades = (($nPositiveColors-1)*($nPositiveShadesPerColor-1));
    
    # negative color information
    $negativeColorString=$negativeColorString.",".$negativeColorString if($negativeColorString !~ /,/);
    my @negativeColorStringArray=split(/,/,$negativeColorString);
    my $nNegativeColors=@negativeColorStringArray;
    my $nNegativeShadesPerColor=floor($nNegativeColorShades/($nNegativeColors-1));
    $nNegativeColorShades = (($nNegativeColors-1)*($nNegativeShadesPerColor-1));
    
    # adjust pos/neg color information to make symmetrical
    $nPositiveColorShades = $nNegativeColorShades if($nNegativeColorShades < $nPositiveColorShades);
    $nNegativeColorShades = $nPositiveColorShades if($nPositiveColorShades < $nNegativeColorShades);
    
    ($colorPalette,$img,$nPositiveColorShades)=calculateColorPalette($nPositiveColorShades,$positiveColorString,\%availableColors,$colorPalette,"pc",$img,$transparency);
    ($colorPalette,$img,$nNegativeColorShades)=calculateColorPalette($nNegativeColorShades,$negativeColorString,\%availableColors,$colorPalette,"nc",$img,$transparency);
        
    my $nPosNegColorShades=$nPositiveColorShades=$nNegativeColorShades;
    my $colorDistance = ($colorScaleEnd-$colorScaleStart);
    my $colorBucketSize = ($colorDistance/($nPosNegColorShades-1));
    
    return($img,$colorPalette,$nPosNegColorShades,$colorBucketSize,\%availableColors);
}

sub drawWebPlot($$$$$$$$$$$$$$) {
    my $matrixObject=shift;
    my $matrix=shift;
    my $minDistance=shift;
    my $maxDistance=shift;
    my $imageWidth=shift;
    my $headerSizing=shift;
    my $colorScaleStartTile=shift;
    my $colorScaleEndTile=shift;
    my $colorScaleStart=shift;
    my $colorScaleEnd=shift;
    my $posColorString=shift;
    my $negColorString=shift;
    my $missingColorString=shift;
    my $drawBGInteractions=shift;
    
    my $inc2header=$matrixObject->{ inc2header };
    my $header2inc=$matrixObject->{ header2inc };
    my $numYHeaders=$matrixObject->{ numYHeaders };
    my $numXHeaders=$matrixObject->{ numXHeaders };
    my $missingValue=$matrixObject->{ missingValue };
    my $output=$matrixObject->{ output };
    my $verbose=$matrixObject->{ verbose };
    
    #globals
    my $includeCis=1;
    my $includeTrans=0;
    my $cisApproximateFactor=1;
    my $excludeZero=0;
    
    # dump matrix data into input lists (CIS + TRANS)
    print STDERR "seperating cis/trans data...\n" if($verbose);
    my ($inputDataCis,$inputDataTrans)=matrix2inputlist($matrixObject,$matrix,$includeCis,$includeTrans,$minDistance,$maxDistance,$excludeZero,$cisApproximateFactor);
    croak "$output - no avaible CIS data" if((scalar @{ $inputDataCis } <= 0) and ($includeCis) and ($includeTrans == 0));
    croak "$output - no avaible TRANS data" if((scalar @{ $inputDataTrans } <= 0) and ($includeTrans) and ($includeCis == 0));
    
    print STDERR "\n" if($verbose);
    
    # calculate the color scale
    if(($colorScaleStart eq "NA") or ($colorScaleEnd eq "NA")) {
        print STDERR "Calculating Auto Color Scale ...\n" if($verbose);
        my ($tmpColorScaleStart,$tmpColorScaleEnd)=autoScale($matrixObject,$matrix,$colorScaleStartTile,$colorScaleEndTile,"all");
        $colorScaleStart=$tmpColorScaleStart if($colorScaleStart eq "NA");
        $colorScaleEnd=$tmpColorScaleEnd if($colorScaleEnd eq "NA");
    }
    print STDERR "\t".$colorScaleStart."-".$colorScaleEnd."\n" if($verbose);
    
    # global image parameters
    my $imageHeight=ceil($imageWidth/3);
    my $numSizes=5;
    my $colorString=$posColorString."___".$negColorString;
    my $binSize=(($colorScaleEnd-$colorScaleStart)/$numSizes);    
    print STDERR "imageHeight\t$imageHeight\n" if($verbose);
    print STDERR "imageWidth\t$imageWidth\n" if($verbose);
    
    print STDERR "\n" if($verbose);
    
    print STDERR "Calculating Color Information...\n" if($verbose);
    my ($img,$colorPalette,$nColorShades,$colorBucketSize)=initImage($colorScaleStart,$colorScaleEnd,$imageHeight,$imageWidth,$colorString);
    print STDERR "\tnColorShades=$nColorShades\n" if($verbose);
    print STDERR "\tcolorBucketSize=$colorBucketSize\n" if($verbose);

    # get left bounds
    my $leftXHeader=$inc2header->{ x }->{0};
    my $leftXHeaderObject=getHeaderObject($leftXHeader);
    my $leftYHeader=$inc2header->{ y }->{0};
    my $leftYHeaderObject=getHeaderObject($leftYHeader);
    
    # get right bounds
    my $rightXHeader=$inc2header->{ x }->{$numXHeaders-1};
    my $rightXHeaderObject=getHeaderObject($rightXHeader);
    my $rightYHeader=$inc2header->{ y }->{$numYHeaders-1};
    my $rightYHeaderObject=getHeaderObject($rightYHeader);
    
    # get region bounds
    my $regionStart = min($leftXHeaderObject->{ start },$leftYHeaderObject->{ start });
    my $regionEnd = max($rightXHeaderObject->{ end },$rightYHeaderObject->{ end });
    my $regionSize=($regionEnd-$regionStart);
    croak "bad region size [$regionSize]" if($regionSize <= 0);
    
    my $heightBorder=floor($imageHeight*.05);
    my $widthBorder=floor($imageWidth*.05);
    my $regionStep=1;
    $regionStep=ceil($regionSize/($imageWidth-($widthBorder*2)));

    # now draw web plot 
    
    # draw the y headers
    for(my $y=0;$y<$numYHeaders;$y++) {
    
        my $yHeader=$inc2header->{ y }->{$y};
        my $yHeaderObject=getHeaderObject($yHeader);
        my $yHeaderStart=$yHeaderObject->{ start };
        my $yHeaderEnd=$yHeaderObject->{ start };
        my $yHeaderName=$yHeaderObject->{ subName };
        my $yHeaderSize=$yHeaderObject->{ size };
        my $yHeaderMidpoint=floor(($yHeaderStart+$yHeaderEnd)/2)-$regionStart;
        my $yHeaderRelativePosition=floor($yHeaderMidpoint/$regionStep);
        
        my $color=$colorPalette->{ G };
        
        # default ball size    (scaled to fragment size)
        my $ballWidth=ceil($yHeaderSize/$regionStep);        
        $ballWidth=10 if($ballWidth > 10);
        $color=$colorPalette->{ G } if(exists($headerSizing->{$yHeader}));
        $ballWidth=20 if(exists($headerSizing->{$yHeader}));
        #next if(!(exists($headerSizing->{$yHeader})));
        
        my @tmp=split(/\|/,$yHeaderName);
        my $name=$tmp[0];
        
        $img->setAntiAliased($color);
        $img->filledArc($yHeaderRelativePosition+$widthBorder,$heightBorder,$ballWidth,$ballWidth,0,360,gdAntiAliased);
        $img->setAntiAliased($colorPalette->{ B });
        $img->arc($yHeaderRelativePosition+$widthBorder,$heightBorder,$ballWidth+1,$ballWidth+1,0,360,gdAntiAliased);
    }
    
    # draw the x headers
    for(my $x=0;$x<$numXHeaders;$x++) {
    
        my $xHeader=$inc2header->{ x }->{$x};
        my $xHeaderObject=getHeaderObject($xHeader);
        my $xHeaderStart=$xHeaderObject->{ start };
        my $xHeaderEnd=$xHeaderObject->{ start };
        my $xHeaderName=$xHeaderObject->{ subName };
        my $xHeaderSize=$xHeaderObject->{ size };
        my $xHeaderMidpoint=floor(($xHeaderStart+$xHeaderEnd)/2)-$regionStart;
        my $xHeaderRelativePosition=floor($xHeaderMidpoint/$regionStep);
        
        my $color=$colorPalette->{ G };
        
        # default ball size    (scaled to fragment size)
        my $ballWidth=ceil($xHeaderSize/$regionStep);        
        $ballWidth=10 if($ballWidth > 10);
        $color=$colorPalette->{ G } if(exists($headerSizing->{$xHeader}));
        $ballWidth=20 if(exists($headerSizing->{$xHeader}));
        #next if(!(exists($headerSizing->{$xHeader})));
    
        my @tmp=split(/\|/,$xHeaderName);
        my $name=$tmp[0];
        
        $img->setAntiAliased($color);
        $img->filledArc($xHeaderRelativePosition+$widthBorder,($imageHeight-$heightBorder),$ballWidth,$ballWidth,0,360,gdAntiAliased);
        $img->setAntiAliased($colorPalette->{ B });
        $img->arc($xHeaderRelativePosition+$widthBorder,($imageHeight-$heightBorder),$ballWidth+1,$ballWidth+1,0,360,gdAntiAliased);
    }

    my $nSkipped=0;
    my $nPlotted=0;
    my $nInt=0;

    $img->setAntiAliased($colorPalette->{ BG });
    $img->setThickness(1);
    
    # draw lines
    for(my $y=0;$y<$numYHeaders;$y++) {
        my $yHeader=$inc2header->{ y }->{$y};
        my $yHeaderObject=getHeaderObject($yHeader);
        my $yHeaderStart=$yHeaderObject->{ start };
        my $yHeaderEnd=$yHeaderObject->{ start };
        my $yHeaderMidpoint=floor(($yHeaderStart+$yHeaderEnd)/2)-$regionStart;
        my $yHeaderRelativePosition=floor($yHeaderMidpoint/$regionStep);
        for(my $x=0;$x<$numXHeaders;$x++) {
            my $xHeader=$inc2header->{ x }->{$x};
            my $xHeaderObject=getHeaderObject($xHeader);
            my $xHeaderStart=$xHeaderObject->{ start };
            my $xHeaderEnd=$xHeaderObject->{ start };            
            my $xHeaderMidpoint=floor(($xHeaderStart+$xHeaderEnd)/2)-$regionStart;
            my $xHeaderRelativePosition=floor($xHeaderMidpoint/$regionStep);
        
            # plot background for all interactions
            $img->line($yHeaderRelativePosition+$widthBorder,$heightBorder,$xHeaderRelativePosition+$widthBorder,($imageHeight-$heightBorder),$colorPalette->{ BG }) if($drawBGInteractions);
        }
    }
    
    # draw lines
    for(my $y=0;$y<$numYHeaders;$y++) {
    
        my $yHeader=$inc2header->{ y }->{$y};
        my $yHeaderObject=getHeaderObject($yHeader);
        my $yHeaderStart=$yHeaderObject->{ start };
        my $yHeaderEnd=$yHeaderObject->{ start };
        my $yHeaderName=$yHeaderObject->{ subName };
        my $yHeaderSize=$yHeaderObject->{ size };
        my $yHeaderMidpoint=floor(($yHeaderStart+$yHeaderEnd)/2)-$regionStart;
        my $yHeaderRelativePosition=floor($yHeaderMidpoint/$regionStep);
        
        for(my $x=0;$x<$numXHeaders;$x++) {
            
            my $xHeader=$inc2header->{ x }->{$x};
            my $xHeaderObject=getHeaderObject($xHeader);
            my $xHeaderStart=$xHeaderObject->{ start };
            my $xHeaderEnd=$xHeaderObject->{ start };
            my $xHeaderName=$xHeaderObject->{ subName };
            my $xHeaderSize=$xHeaderObject->{ size };
            my $xHeaderMidpoint=floor(($xHeaderStart+$xHeaderEnd)/2)-$regionStart;
            my $xHeaderRelativePosition=floor($xHeaderMidpoint/$regionStep);
            
            $nInt++;
            
            if(exists($matrix->{$y}->{$x})) {
                my $score=$matrix->{$y}->{$x};
            
                next if($score eq "NA");
                
                my $sizeIndex = getColorIndex($score,$colorScaleStart,$colorScaleEnd,$numSizes,$binSize);
                my $colorIndex = getColorIndex($score,$colorScaleStart,$colorScaleEnd,$nColorShades,$colorBucketSize);
                
                # hide all signal below color scale start
                next if($colorIndex == 0);
                
                $img->setThickness($sizeIndex+1);
                
                my $color="NULL";
                if($score eq "NA") {
                    $color=$colorPalette->{$missingColorString};
                } elsif($score < 0) {
                    croak "colorIndex [$colorIndex] [nc] does not exists! ($score - $colorIndex)" if(!(exists($colorPalette->{ nc }[$colorIndex])));
                    $color=$colorPalette->{ nc }[$colorIndex];
                } elsif($score >= 0) {
                    croak "colorIndex [$colorIndex] [pc] does not exists! ($score - $colorIndex)" if(!(exists($colorPalette->{ nc }[$colorIndex])));
                    $color=$colorPalette->{ pc }[$colorIndex];
                }
                
                $nPlotted++;
                
                $img->line($yHeaderRelativePosition+$widthBorder,$heightBorder,$xHeaderRelativePosition+$widthBorder,($imageHeight-$heightBorder),$color);
            }
            
        }
    }

    open(OUT,">".$output.".webplot.png") or croak "Could not open file [$output] - $!";
    binmode OUT;
    print OUT $img->png;
    close(OUT);

}

sub getHeaderSizing($) {
    my $headerSizingFile=shift;

    croak "headerSizingFile [$headerSizingFile] does not exist" if(!(-e $headerSizingFile));

    my %headerSizing=();
    
    open(IN,inputWrapper($headerSizingFile)) or croak "Could not open file [$headerSizingFile] - $!";
    while(my $line = <IN>) {
        chomp($line);
        next if(($line eq "") or ($line =~ m/^#/));
        
        my @tmp=split(/\t/,$line);
        my $name=$tmp[3];
        my $score=1;
        $score=$tmp[4] if(defined($tmp[4]));
        
        next if($score == 0);
        
        $headerSizing{$name}=$score;
        
    }
    
    return(\%headerSizing);
}

my %options;
my $results = GetOptions( \%options,'inputMatrix|i=s','verbose|v','output|o=s','minDistance|minDist=i','maxDistance|maxDist=i','imageWidth|iw=s','headerSizingFile|hsf=s','colorScaleStartTile|startTile=s','colorScaleEndTile|endTile=s','colorScaleStart|start=s','colorScaleEnd|end=s','posColorString|pc=s','negColorString|nc=s','missingColorString|mc=s','drawBGInteractions|dbi') or croak help();
my ($ret,$inputMatrix,$verbose,$output,$minDistance,$maxDistance,$imageWidth,$headerSizingFile,$colorScaleStartTile,$colorScaleEndTile,$colorScaleStart,$colorScaleEnd,$posColorString,$negColorString,$missingColorString,$drawBGInteractions)=check_options( \%options );

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
my $numHeaders=$matrixObject->{ numTotalHeaders };
my $missingValue=$matrixObject->{ missingValue };
my $contigList=$matrixObject->{ contigList };
my $numContigs=$matrixObject->{ numContigs };
my $index2contig=$matrixObject->{ index2contig };
my $symmetrical=$matrixObject->{ symmetrical };
my $equalHeaderFlag=$matrixObject->{ equalHeaderFlag };
my $inputMatrixName=$matrixObject->{ inputMatrixName };
$output=$matrixObject->{ output };

croak "matrix must be symmetrical" if($symmetrical == 0);

print STDERR "contigs:\n" if($verbose);
for(my $c=0;$c<$numContigs;$c++) {
    my $contig=$index2contig->{ xy }->{$c};
    my $contigStart=$contigList->{$contig}->{ contigStart };
    my $contigEnd=$contigList->{$contig}->{ contigEnd };
    my $contigAssembly=$contigList->{$contig}->{ contigAssembly };
    my $contigChromosome=$contigList->{$contig}->{ contigChromosome };
    my $contigLength=$contigList->{$contig}->{ contigLength };
    print STDERR "\t$contig\t$contigAssembly\t$contigChromosome\t$contigStart - $contigEnd [$contigLength]\n" if($verbose);
}
    
print STDERR "\n" if($verbose);

#read Matrix
my ($matrix)=getData($inputMatrix,$matrixObject,$verbose);

print STDERR "\n" if($verbose);

my $headerSizing={};
$headerSizing=getHeaderSizing($headerSizingFile) if($headerSizingFile ne "");

drawWebPlot($matrixObject,$matrix,$minDistance,$maxDistance,$imageWidth,$headerSizing,$colorScaleStartTile,$colorScaleEndTile,$colorScaleStart,$colorScaleEnd,$posColorString,$negColorString,$missingColorString,$drawBGInteractions);

print STDERR "\n" if($verbose);
