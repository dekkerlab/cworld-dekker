<img height=40 src='http://my5C.umassmed.edu/images/3DG.png' title='3D-Genome' />
&nbsp;&nbsp;
<img height=30 src='http://my5C.umassmed.edu/images/dekkerlabbioinformatics.gif' />
&nbsp;&nbsp;
<img height=30 src='http://my5C.umassmed.edu/images/umasslogo.gif' />

# cworld::dekker

a collection of perl/python/R scripts for maniputing 3C/4C/5C/Hi-C data.

## Installation

This package requires:
```
libgd (2.0.28 or higher)
bedtools (bedtools or bedtools2)
R / Rscript
```

To download/install dependencies:
```
libgd - https://github.com/libgd/libgd/releases
bedtools - https://github.com/arq5x/bedtools2/releases
R - https://cran.r-project.org/bin/

Please follow installation instructions per dependency.
```

LIBGD Troubleshooting 

```
For Mac OSX:
brew install gd

For Linux:
yum install libgd
```

Download the project.
```
wget https://github.com/blajoie/cworld-dekker/archive/master.zip
```
or clone the git project
```
[ssh] - git clone git@github.com:blajoie/cworld-dekker.git
[https] - git clone https://github.com/blajoie/cworld-dekker.git
```

Unzip the master:
```
unzip cworld-dekker-master.zip
cd cworld-dekker-master/
```

To install the module:
```
perl Build.PL
./Build
./Build install
```

After installing the module, you should be free to run the any scripts found in scriots/
e.g.
```
$ perl scripts/heatmap.pl
```

## Full Documentation

WIP

## Communication

- [Bryan Lajoie](https://github.com/blajoie)
- Twitter: [@my5C](https://twitter.com/my5C)

## What does it do?

cworld::dekker is a large collection of perl / python / R scripts for manipulating 3C/4C/5C/Hi-C data

## Tools

```
Tool:		addMatrixHeaders.pl
Version:	0.01
Summary:	add headers to a matrix txt file

Usage: perl addMatrixHeaders.pl [OPTIONS] -i <inputMatrix> -xhf <xHeaderFile> -yhf <yHeaderFile>

Required:
	-i         []         input matrix file
	--xhf      []         x-axis header file
	--yhf      []         y-axis header file

Options:
	-v         []         FLAG, verbose mode
	-o         []         prefix for output file(s)

Notes:
    This script can add headers to a _naked_ matrix file.
    Matrix can be TXT or gzipped TXT.
    See website for matrix format details.

Contact:
    Bryan R. Lajoie
    Dekker Lab 2015
    https://github.com/blajoie/cworld-dekker
    http://my5C.umassmed.edu
```


```
Tool:		aggregrateBED.pl
Version:	0.01
Summary:	sliding window aggregrate of BED5 data

Usage: perl aggregrateBED.pl [OPTIONS] -i <inputBEDFile> -a <assembly>

Required:
	-i         []         input BED5 file
	--a        []         assembly of BED5 file

Options:
	-v         []         FLAG, verbose mode
	-o         []         prefix for output file(s)
	--cbf      []         contig bound file (bed3+), to specify start/end of regions/contigs.
	--wsize    [100000]   window size, in bp
	--wstep    [1]        window step factor. 1 = non overlapping windows, 2-inf = overlapping windows
	--wmode    [sum]      window signal aggregrate mode.  available modes are sum,mean,median,binary,min,max,stdev,variance etc
	--ms       [NA]       suplment missing scores with supplied value. (NA,0)
	--cis      []         enable cis mode - seperate track for each contig
	--ez       []         FLAG, ignore 0s in all calculations
	--yl       [-1]       optional ylimit for plot, by default, plots are autoScaled

Notes:
    This program can window a BED5 formatted file into fixed size intervals.
    BED signal is then aggregrated (mean,median,sum etc)

Contact:
    Bryan R. Lajoie
    Dekker Lab 2015
    https://github.com/blajoie/cworld-dekker
    http://my5C.umassmed.edu
```


```
Tool:		anchorPileUp.pl
Version:	0.01
Summary:	pile up cData around specified list of 'anchors'

Usage: perl anchorPileUp.pl [OPTIONS] -i <inputMatrix> -abf <anchorBedFile>

Required:
	-i         []         input matrix file
	--abf      []         anchor bed file

Options:
	-v         []         FLAG, verbose mode
	-o         []         prefix for output file(s)
	--lof      []         optional loess object file (pre-calculated loess)
	--azs      []         anchorZoneSize, size of zone to use around anchor (x axis - in BP)
	--am       []         aggregrate mode, how to aggregrate signal [mean,sum,median]
	--minDist  []         minimum allowed interaction distance, exclude < N distance (in BP)
	--maxDist  []         maximum allowed interaction distance, exclude > N distance (in BP)
	--ez       []         FLAG, ignore 0s in all calculations

Notes:
    This script aggregrates cData surrounding a list of 'anchors'.
    Input Matrix can be TXT or gzipped TXT.
    AnchorListFile must be BED5+ format.
    See website for matrix format details.

Contact:
    Bryan R. Lajoie
    Dekker Lab 2015
    https://github.com/blajoie/cworld-dekker
    http://my5C.umassmed.edu
```


```
Tool:		anchorPurge.pl
Version:	0.01
Summary:	filters out row/col from C data matrix file

Usage: perl anchorPurge.pl [OPTIONS] -i <inputMatrix>

Required:
	-i         []         input matrix file

Options:
	-v         []         FLAG, verbose mode
	-o         []         prefix for output file(s)
	--lof      []         optional loess object file (pre-calculated loess)
	--ta       [0.9]      fraction of data to keep, 0.9 = trim 10% off
	--ic       []         model CIS data to detect outlier row/cols
	--it       []         model TRANS data to detect outlier row/cols
	--mof      []         optional manual outlier file, 1 header per line to be filtered.
	--minDist  []         minimum allowed interaction distance, exclude < N distance (in BP)
	--maxDist  []         maximum allowed interaction distance, exclude > N distance (in BP)
	--ca       [0.01]     lowess alpha value, fraction of datapoints to smooth over
	--dif      []         FLAG, disable loess IQR (outlier) filter
	--caf      [1000]     cis approximate factor to speed up loess, genomic distance / -caffs
	--lt       []         log transform input data into specified base
	--ez       []         FLAG, ignore 0s in all calculations
	--fm       []         outlier detection mode - zScore,obsExp

Notes:
    This script detects and removes outlier row/cols.
    Outlier row/cols can be detected by either cis or trans data.
    Matrix can be TXT or gzipped TXT.
    See website for matrix format details.

Contact:
    Bryan R. Lajoie
    Dekker Lab 2015
    https://github.com/blajoie/cworld-dekker
    http://my5C.umassmed.edu
```


```
Tool:		applyCorrection.pl
Version:	0.01
Summary:	apply correction to factor - external factors

Usage: perl applyCorrection.pl [OPTIONS] -i <inputMatrix>

Required:
	-i         []         input matrix file
	--ff       []         row/col factor file

Options:
	-v         []         FLAG, verbose mode
	-o         []         prefix for output file(s)
	--lof      []         optional loess object file (pre-calculated loess)
	--ic       []         model CIS data to detect outlier row/cols
	--it       []         model TRANS data to detect outlier row/cols
	--minDist  []         minimum allowed interaction distance, exclude < N distance (in BP)
	--maxDist  []         maximum allowed interaction distance, exclude > N distance (in BP)
	--ca       [0.01]     lowess alpha value, fraction of datapoints to smooth over
	--dif      []         FLAG, disable loess IQR (outlier) filter
	--caf      [1000]     cis approximate factor to speed up loess, genomic distance / -caffs
	--lt       []         log transform input data into specified base
	--ez       []         FLAG, ignore 0s in all calculations
	--fm       []         outlier detection mode - zScore,obsExp

Notes:
    This script can apply a set of primer factors to a matrix [normalize].
    Matrix can be TXT or gzipped TXT.
    See website for matrix format details.

Contact:
    Bryan R. Lajoie
    Dekker Lab 2015
    https://github.com/blajoie/cworld-dekker
    http://my5C.umassmed.edu
```


```
Tool:		binMatrix.pl
Version:	0.01
Summary:	bin/aggregrate matrix into fixed sized intervals

Usage: perl binMatrix.pl [OPTIONS] -i <inputMatrix> -bsize <binSize> -bstep <binStep> -bmode <binMode>

Required:
	-i         []         input matrix file
	-o         []         prefix for output file(s)
	--bsize    []         bin size in bp
	--bstep    []         bin step in factor of overlap, 1=non overlapping, 2-inf=overlapping
	--bmode    []         bin mode, how to aggregrate data, mean,sum,median etc.

Options:
	-v         []         FLAG, verbose mode
	--ez       []         FLAG, ignore 0s in all calculations
	--nc       []         FLAG, use nice(clean) coordinates.  start first bin at position cleanly divisible by window size parameter.

Notes:
    This script can bin a matrix into fixed size intervals and aggregrate the data.
    Matrix can be TXT or gzipped TXT.
    See website for matrix format details.

Contact:
    Bryan R. Lajoie
    Dekker Lab 2015
    https://github.com/blajoie/cworld-dekker
    http://my5C.umassmed.edu
```


```
Tool:		changeMatrixHeaders.pl
Version:	0.01
Summary:	replace matrix row/col headers, or subset matrix by list of headers

Usage: perl changeMatrixHeaders.pl [OPTIONS] -i <inputMatrix> -hmf <headerMapFile>

Required:
	-i         []         input matrix file
	--hmf      []         header map file, old (tab) new

Options:
	-v         []         FLAG, verbose mode
	-o         []         prefix for output file(s)

Notes:
    This script can bin a matrix into fixed size intervals and aggregrate the data.
    Matrix can be TXT or gzipped TXT.
    See website for matrix format details.

Contact:
    Bryan R. Lajoie
    Dekker Lab 2015
    https://github.com/blajoie/cworld-dekker
    http://my5C.umassmed.edu
```


```
Tool:		collapseMatrix.pl
Version:	0.01
Summary:	collapse matrix by (chr,name,group), sum signal

Usage: perl collapseMatrix.pl [OPTIONS] -i <inputMatrix>

Required:
	-i         []         input matrix file

Options:
	-v         []         FLAG, verbose mode
	-o         []         prefix for output file(s)
	--cb       [chr]      method to collapse headers (chr,name,group)
	--ed       []         FLAG, exclude diagonal

Notes:
    This script collapses headers by chr,name,group.
    Matrix can be TXT or gzipped TXT.
    See website for matrix format details.

Contact:
    Bryan R. Lajoie
    Dekker Lab 2015
    https://github.com/blajoie/cworld-dekker
    http://my5C.umassmed.edu
```


```
Tool:		column2matrix.pl
Version:	0.01
Summary:	turn list (3 tab) file into matrix

Usage: perl column2matrix.pl [OPTIONS] -i <inputMatrix>

Required:
	-i         []         input list file

Options:
	-v         []         FLAG, verbose mode
	-o         []         prefix for output file(s)
	--oxh      []         optional x-axis header file
	--oyh      []         optional y-axis header file
	--mv       []         optional missing-value [0,NA]

Notes:
    This script converts a 3 column txt (tsv) file to a matrix.
    Use -oxh and -oyh to set the x-axis / y-axis headers.
    See website for matrix format details.

Contact:
    Bryan R. Lajoie
    Dekker Lab 2015
    https://github.com/blajoie/cworld-dekker
    http://my5C.umassmed.edu
```


```
Tool:		combineMatrices.pl
Version:	0.01
Summary:	combine matrices [sum,mean,median,min,max]

Usage: perl combineMatrices.pl [OPTIONS] -i <inputMatrix>

Required:
	-i@        []         input matrix file array [MULTIPLE]

Options:
	-v         []         FLAG, verbose mode
	-o         []         prefix for output file(s)
	--n        []         optional, prefix for output file
	--cm       []         optional, combine mode [sum,mean,median,min,max]
	--tmp      [/tmp]     optional, tmp direction for tmp files

Notes:
    This script can combine multiple matrices.
    See website for matrix format details.

Contact:
    Bryan R. Lajoie
    Dekker Lab 2015
    https://github.com/blajoie/cworld-dekker
    http://my5C.umassmed.edu
```


```
Tool:		compareInsulation.pl
Version:	0.01
Summary:	compare insulation vector - calculate difference

Usage: perl compareInsulation.pl [OPTIONS] -i1 <insulation_1> -i2 <insulation_2>

Required:
	-1         []         insulation vectror 1 file
	-2         []         insulation vectror 2 file

Options:
	-v         []         FLAG, verbose mode
	-o         []         prefix for output file(s)
	--yb       []         optional, fix y-limit for insulation image
	--bg       []         optional, enable transparent background for insulation image

Notes:
    This script can compare two insulation vectors [difference].
    See website for matrix format details.

Contact:
    Bryan R. Lajoie
    Dekker Lab 2015
    https://github.com/blajoie/cworld-dekker
    http://my5C.umassmed.edu
```


```
Tool:		compareMatrices.pl
Version:	0.01
Summary:	performs comparison between two matrices

Usage: perl compareMatrices.pl [OPTIONS] -i1 <inputMatrix_1> -i2 <inputMatrix_2>

Required:
	-1         []         input matrix 1 file
	-2         []         input matrix 2 file

Options:
	-v         []         FLAG, verbose mode
	-o         []         prefix for output file(s)
	--cm       []         optional, compare mode [log2ratio,add,sum,mean,subtract,divide,multiply,min,max]

Notes:
    This script can compare two matrices.
    See website for matrix format details.

Contact:
    Bryan R. Lajoie
    Dekker Lab 2015
    https://github.com/blajoie/cworld-dekker
    http://my5C.umassmed.edu
```


```
Tool:		correlateMatrices.pl
Version:	0.01
Summary:	performs correlation between two matrices

Usage: perl correlateMatrices.pl [OPTIONS] -i1 <inputMatrix_1> -i2 <inputMatrix_2>

Required:
	-1         []         input matrix 1 file
	-2         []         input matrix 2 file

Options:
	-v         []         FLAG, verbose mode
	-o         []         prefix for output file(s)
	--ppd      []         FLAG, plot correlation per distance
	--ec       []         FLAG, exclude CIS data
	--et       []         FLAG, exclude TRANS data
	--ez       []         FLAG, exclude zeros from all calculations
	--ed       []         FLAG, exclude diagonal bin (y=x) from all calculations
	--minDist  []         minimum allowed interaction distance, exclude < N distance (in BP)
	--maxDist  []         maximum allowed interaction distance, exclude > N distance (in BP)
	--cm       []         optional, correlation mode (pearson,spearman)
	--lt       []         optional, log transform data (2,10)
	--fbp      []         optional, fixed bin population size (num of interactions) (-1=auto, 0= ff)
	--of       []         optional, outlierFraction [0-1], removed from top/bottom during correlation
	--ymin     []         optional, y min value of plot
	--ymax     []         optional, y max value of plot
	--xmin     []         optional, x min value of plot
	--xmax     []         optional, x max value of plot
	--tmp      [/tmp]     optional, tmp direction for tmp files

Notes:
    This script can compare two matrices via scatter plot regression.
    See website for matrix format details.

Contact:
    Bryan R. Lajoie
    Dekker Lab 2015
    https://github.com/blajoie/cworld-dekker
    http://my5C.umassmed.edu
```


```
Tool:		coverageCorrect.pl
Version:	0.01
Summary:	can perform coverage correction on matrix [balancing] cis/trans

Usage: perl coverageCorrect.pl [OPTIONS] -i <inputMatrix>

Required:
	-i         []         input matrix file

Options:
	-v         []         FLAG, verbose mode
	-o         []         prefix for output file(s)
	--ic       []         model CIS data to detect outlier row/cols
	--it       []         model TRANS data to detect outlier row/cols
	--minDist  []         minimum allowed interaction distance, exclude < N distance (in BP)
	--maxDist  []         maximum allowed interaction distance, exclude > N distance (in BP)
	--ca       [0.01]     lowess alpha value, fraction of datapoints to smooth over
	--dif      []         FLAG, disable loess IQR (outlier) filter
	--caf      [1000]     cis approximate factor to speed up loess, genomic distance / -caffs
	--lt       []         log transform input data into specified base
	--ez       []         FLAG, ignore 0s in all calculations
	--fm       [zScore]   outlier detection mode - zScore,obsExp
	--ct       [0.1]      optional, convergance threshold

Notes:
    This script can compare two matrices via scatter plot regression.
    See website for matrix format details.

Contact:
    Bryan R. Lajoie
    Dekker Lab 2015
    https://github.com/blajoie/cworld-dekker
    http://my5C.umassmed.edu
```


```
Tool:		digitizePicture.pl
Version:	0.01
Summary:	digitize picture into my5C matrix format

Usage: perl digitizePicture.pl [OPTIONS] -i <inputPNG>

Required:
	-i         []         input matrix file

Options:
	-v         []         FLAG, verbose mode
	-o         []         prefix for output file(s)
	--cm       [mean]     color mode (red,blue,green,min,max,mean,median)

Notes:
    This script turns any PNG into a my5C formatted matrix.

Contact:
    Bryan R. Lajoie
    Dekker Lab 2015
    https://github.com/blajoie/cworld-dekker
    http://my5C.umassmed.edu
```


```
Tool:		extractSubMatrices.pl
Version:	0.01
Summary:	this script can extract sub matrices classified by chr/group/name

Usage: perl extractSubMatrices.pl [OPTIONS] -i <inputMatrix>

Required:
	-i         []         input matrix file

Options:
	-v         []         FLAG, verbose mode
	-o         []         prefix for output file(s)
	--ozc      []         @ [MULTIPLE], genomic coordinates for subset (-ozc chr:start-end)
	--eb       [chr]      optional, classifier by which to extract sub matrices (chr,region,group)
	--eco      []         FLAG, extract only cis matrices
	--usn      []         FLAG, use short input file names
	--tmp      [/tmp]     optional, tmp direction for tmp files

Notes:
    This script extracts all sub matrices of a given input matrix.
    Matrix can be TXT or gzipped TXT.
    See website for matrix format details.

Contact:
    Bryan R. Lajoie
    Dekker Lab 2015
    https://github.com/blajoie/cworld-dekker
    http://my5C.umassmed.edu
```


```
Tool:		fillMissingData.pl
Version:	0.01
Summary:	replace NAs with expected signals

Usage: perl fillMissingData.pl [OPTIONS] -i <inputMatrix>

Required:
	-i         []         input matrix file

Options:
	-v         []         FLAG, verbose mode
	-o         []         prefix for output file(s)
	--lof      []         optional loess object file (pre-calculated loess)
	--ec       []         FLAG, exclude CIS data
	--et       []         FLAG, exclude TRANS data
	--ca       [0.01]     lowess alpha value, fraction of datapoints to smooth over
	--dif      []         FLAG, disable loess IQR (outlier) filter
	--caf      [1000]     cis approximate factor to speed up loess, genomic distance / -caffs
	--minDist  []         minimum allowed interaction distance, exclude < N distance (in BP)
	--maxDist  []         maximum allowed interaction distance, exclude > N distance (in BP)
	--ez       []         FLAG, ignore 0s in all calculations

Notes:
    This script can fill NAs with the lowess expected.
    Matrix can be TXT or gzipped TXT.
    See website for matrix format details.

Contact:
    Bryan R. Lajoie
    Dekker Lab 2015
    https://github.com/blajoie/cworld-dekker
    http://my5C.umassmed.edu
```


```
Tool:		generateBins.pl
Version:	0.01
Summary:	create my5C formatted headers

Usage: perl generateBins.pl [OPTIONS] -a <assembly> -r <regionCoordinates -bsize <binSize> -bstep <binStep>

Required:
	--a        [1]        genome assembly
	--r        []         region coordinates, e.g. chr1:1-50000000 (UCSC)
	--bsize    []         bin size in bp
	--bstep    []         bin step in factor of overlap, 1=non overlapping, 2-inf=overlapping

Options:
	-v         []         FLAG, verbose mode
	--n        [myHeaders] optional, prefix for output files

Notes:
    This script can create my5C formatted headers.
    Matrix can be TXT or gzipped TXT.
    See website for matrix format details.

Contact:
    Bryan R. Lajoie
    Dekker Lab 2015
    https://github.com/blajoie/cworld-dekker
    http://my5C.umassmed.edu
```


```
Tool:		heatmap.pl
Version:	0.01
Summary:	draws heatmap PNG of matrix file

Usage: perl heatmap.pl [OPTIONS] -i <inputMatrix>

Required:
	-i@        []         input matrix file [MULTIPLE]

Options:
	-v         []         FLAG, verbose mode
	-o         []         prefix for output file(s)
	--sfs      []         FLAG, scaleFragmentSizes, scale pixels by fragment/bin size
	--dt       []         FLAG, drawTriangle, draw heatmap as triangle (only upper triangle)
	--dd       []         FLAG, drawDiamond, draw heatmap as diamond
	--bg       []         FLAG, use transparent background
	--em       []         FLAG, embed meta data into PNG (warning - usage on large images (>5000pixels) will use ~GB mem
	--dpb      []         FLAG, draw pixel border
	--ocb      []         FLAG, omit contig/region border
	--maxDist  []         FLAG, draw header labels
	--ds       []         FLAG, draw scores in pixels
	--name     []         prefix for output file
	--is       [800]      ideal image size in pixel (larger of width/height
	--ps       []         x/y axis pixel size
	--yps      []         y axis pixel size
	--xps      []         x axis pixel size
	--maxdim   [30000]    maximum image dimension [hard-limited:32000]
	--lt       [0]        log transform input data into specified base
	--start    []         absolute value for color start
	--end      []         absolute vlaue for color end
	--startTile [0.025]    fraction value for color start
	--endTile  [0.975]    fraction value for color end
	--hbf      []         highlight bed file - overlap row/col will be highlighted
	--iq       [9]        image quality, 0=best 9=worst
	--sm       []         scale mode, single color scale, or seperate for cis/trans [combined/seperate]
	--pc       [white,orange,red,darkRed] positive color scale string
	--nc       [white,cyan,blue,darkBlue] negative color scale string
	--mc       [null]     missing data color
	--hc       [cyan]     highlight row/col data color
	--cs       []         contig spacing, pixel size for spacing between contigs (contig border)

Notes:
    This script coverts a matrix into a PNG heatmap image.
    Matrix should have row/col headers in the standard my5C format.
    Matrix can be TXT or gzipped TXT.
    See website for matrix format details.

Contact:
    Dekker Lab
    http://my5C.umassmed.edu
    my5C.help@umassmed.edu
```


```
Tool:		insulation2tads.pl
Version:	0.01
Summary:	create tad specific headers

Usage: perl insulation2tads.pl [OPTIONS] -i <inputMatrix> -b <boundaryFile>

Required:
	-i         []         cWorld insulation file
	--b        []         cWorld boundary file

Options:
	-v         []         FLAG, verbose mode
	-o         []         prefix for output file(s)
	--mbs      []         min boundary strength (0-inf)
	--mts      []         min tad strength (0-inf)
	--yb       [auto]     yBound, 0 - yBound for bedGraph

Notes:
    This script can assemble consecutive non-overlapping TADs from a insulation/boundary file.

Contact:
    Bryan R. Lajoie
    Dekker Lab 2015
    https://github.com/blajoie/cworld-dekker
    http://my5C.umassmed.edu
```


```
Tool:		interactionPileUp.pl
Version:	0.01
Summary:	pile up cData around specified list of 'elements'

Usage: perl interactionPileUp.pl [OPTIONS] -i <inputMatrix> -ebf <elementBedFile>

Required:
	-i         []         input matrix file
	--ebf@     []         element bed file (can accept multiple files)

Options:
	-v         []         FLAG, verbose mode
	-o         []         prefix for output file(s)
	--en       []         elementName, descriptor for output files
	-zs        [100000]   elementZoneSize, size of zone to use around element (x axis - in BP)
	--am       []         aggregrateMode, how to aggregrate data (mean,sum,median,iqrMean,min,max etc)
	--minDist  []         minimum allowed interaction distance, exclude < N distance (in BP)
	--maxDist  []         maximum allowed interaction distance, exclude > N distance (in BP)
	--minED    []         minElementDistance, minimum element distance (#elements) to consider
	--maxED    []         maxElementDistance, maximum element distance (#elements) to consider
	--ez       []         FLAG, exclude 0s in all calculations
	--id       [off]      FLAG, include any interaction space that touches diagonal
	--ec       []         FLAG, exclude CIS interactions
	--et       []         FLAG, exclude TRANS interactions
	--d        []         FLAG, debug mode (print out all element-element matrices)

Notes:
    This script aggregrates cData surrounding a list of 'elements'.
    Input Matrix can be TXT or gzipped TXT.
    AnchorListFile must be BED5+ format.
    See website for matrix format details.

Contact:
    Dekker Lab
    http://my5C.umassmed.edu
    my5C.help@umassmed.edu
```


```
Tool:		matrix2anchorPlot.pl
Version:	0.01
Summary:	transform each row/col into 4C style 'anchor' plot.

Usage: perl matrix2anchorPlot.pl [OPTIONS] -i <inputMatrix>

Required:
	-i         []         input matrix file

Options:
	-v         []         FLAG, verbose mode
	-o         []         prefix for output file(s)
	--lof      []         optional, pre-calculated *.loess file
	--ca       [0.01]     lowess alpha value, fraction of datapoints to smooth over
	--dif      []         FLAG, disable loess IQR (outlier) filter
	--caf      [1000]     cis approximate factor to speed up loess, genomic distance / -caffs
	--ez       []         FLAG, ignore 0s in all calculations
	--minDist  []         minimum allowed interaction distance, exclude < N distance (in BP)
	--maxDist  []         maximum allowed interaction distance, exclude > N distance (in BP)
	--hm       []         highlight matrix file - binary flag for specific interactions
	--ymin     []         optional, ymin value for anchor plot y-axis
	--ymax     []         optional, ymax value for anchor plot y-axis
	--slf      []         optional, list of headers to draw, all other skipped

Notes:
    This script can transform each row/col into 4C style 'anchor' plot.
    Matrix can be TXT or gzipped TXT.
    See website for matrix format details.

Contact:
    Bryan R. Lajoie
    Dekker Lab 2015
    https://github.com/blajoie/cworld-dekker
    http://my5C.umassmed.edu
```


```
Tool:		matrix2bed12.pl
Version:	0.01
Summary:	transform matrix into bed12 format (track per row)

Usage: perl matrix2bed12.pl [OPTIONS] -i <inputMatrix>

Required:
	-i         []         input matrix file

Options:
	-v         []         FLAG, verbose mode
	-o         []         prefix for output file(s)
	--start    []         absolute value for display start
	--end      []         absolute vlaue for display end
	--startTile [0.025]    fraction value for display start
	--endTile  [0.975]    fraction value for display end

Notes:
    This script transform a matrix into bed12 format (1 track per row/col).
    Matrix can be TXT or gzipped TXT.
    See website for matrix format details.

Contact:
    Bryan R. Lajoie
    Dekker Lab 2015
    https://github.com/blajoie/cworld-dekker
    http://my5C.umassmed.edu
```


```
Tool:		matrix2compartment.pl
Version:	0.01
Summary:	perform PCA on input matrix

Usage: perl matrix2compartment.pl [OPTIONS] -i <inputMatrix>

Required:
	-i         []         input matrix file

Options:
	-v         []         FLAG, verbose mode
	-o         []         prefix for output file(s)
	--lof      []         optional loess object file (pre-calculated loess)
	--ec       []         FLAG, exclude CIS data
	--et       []         FLAG, exclude TRANS data
	--ca       [0.01]     lowess alpha value, fraction of datapoints to smooth over
	--caf      [1000]     cis approximate factor to speed up loess, genomic distance / -caffs
	--dif      []         FLAG, disable loess IQR (outlier) filter
	--minDist  []         minimum allowed interaction distance, exclude < N distance (in BP)
	--maxDist  []         maximum allowed interaction distance, exclude > N distance (in BP)
	--ez       []         FLAG, ignore 0s in all calculations

Notes:
    This script performs a PCA analysis on the input matrix.
    Eigen1,Eigen2 and Eigen3 are plotted.
    The positive Eigen1 'group' is defined as the most gene rich 'group'.
    AnchorListFile must be BED5+ format.
    See website for matrix format details.

Contact:
    Bryan R. Lajoie
    Dekker Lab 2015
    https://github.com/blajoie/cworld-dekker
    http://my5C.umassmed.edu
```


```
Tool:		matrix2direction.pl
Version:	0.01
Summary:	calculate directionality [tads] on matrix

Usage: perl matrix2direction.pl [OPTIONS] -i <inputMatrix>

Required:
	-i         []         input matrix file

Options:
	-v         []         FLAG, verbose mode
	-o         []         prefix for output file(s)
	--dm       [mean]     optional, directionality aggregrate mode [mean,median]
	--ds       [500000]   optional, directionality window size in bp
	--ez       []         FLAG, ignore 0s in all calculations
	--yb       [auto]     optional, -yBound - +yBound for direction plot y-axis range
	--bg       []         FLAG, use transparent background of direction plot

Notes:
    This script can calculate the directionality log2[(eft/right) for each locus.
    Matrix can be TXT or gzipped TXT.
    See website for matrix format details.

Contact:
    Bryan R. Lajoie
    Dekker Lab 2015
    https://github.com/blajoie/cworld-dekker
    http://my5C.umassmed.edu
```


```
Tool:		matrix2distance.pl
Version:	0.01
Summary:	cumlative reads versus distance

Usage: perl matrix2distance.pl [OPTIONS] -i <inputMatrix>

Required:
	-i         []         input matrix file

Options:
	-v         []         FLAG, verbose mode
	-o         []         prefix for output file(s)
	--minDist  []         minimum allowed interaction distance, exclude < N distance (in BP)
	--maxDist  []         maximum allowed interaction distance, exclude > N distance (in BP)
	--lt       []         optional, log transform data (2,10)
	--ez       []         FLAG, exclude zeros from all calculations
	--ed       []         FLAG, exclude diagonal bin (y=x) from all calculations
	--sn       []         FLAG, skip NAs in calculation

Notes:
    This script can compute the cumlative reads across distances for a given matrix.
    See website for matrix format details.

Contact:
    Bryan R. Lajoie
    Dekker Lab 2015
    https://github.com/blajoie/cworld-dekker
    http://my5C.umassmed.edu
```


```
Tool:		matrix2headerBed.pl
Version:	0.01
Summary:	dump matrix headers as BED file

Usage: perl matrix2headerBed.pl [OPTIONS] -i <inputMatrix>

Required:
	-i         []         input matrix file

Options:
	-v         []         FLAG, verbose mode
	-o         []         prefix for output file(s)

Notes:
    This script dumps a my5C formatted matrix headers into a BED5 file.
    Input Matrix can be TXT or gzipped TXT.
    See website for matrix format details.

Contact:
    Bryan R. Lajoie
    Dekker Lab 2015
    https://github.com/blajoie/cworld-dekker
    http://my5C.umassmed.edu
```


```
Tool:		getMatrixInfo.pl
Version:	0.01
Summary:	get matrix info

Usage: perl getMatrixInfo.pl [OPTIONS] -i <inputMatrix>

Required:
	-i         []         input matrix file

Options:
	-v         []         FLAG, verbose mode
	-o         []         prefix for output file(s)
	--id       []         FLAG, ignore diagonal bin during normalization

Notes:
    This script returns matrix info (# row/cols, # reads, cis/trans etc).
    Matrix should have row/col headers in the standard my5C format.
    Matrix can be TXT or gzipped TXT.
    See website for matrix format details.

Contact:
    Bryan R. Lajoie
    Dekker Lab 2015
    https://github.com/blajoie/cworld-dekker
    http://my5C.umassmed.edu
```


```
Tool:		matrix2insulation.pl
Version:	0.01
Summary:	calculate insulation index (TADs) of supplied matrix

Usage: perl matrix2insulation.pl [OPTIONS] -i <inputMatrix>

Required:
	-i         []         input matrix file

Options:
	-v         []         FLAG, verbose mode
	-o         []         prefix for output file(s)
	--is       []         optional, insulation square size, size of the insulation square in bp
	--ss       []         optional, insulation smooth size, size of insulation vector smoothing vector in bp
	--ass      []         optional, insulation amplitude spread size, size of amplitude stdev smoothing vector in bp
	--ids      []         optional, insulation delta span (window size of insulationd delta - blue line)
	--im       []         optional, insulation mode (how to aggregrate signal within insulation square), mean,sum,median
	--nt       [0.1]      optional, noise threshold, minimum depth of valley
	--bmoe     [3]        optional, boundary margin of error (specified in number of BINS), added to each side of the boundary
	--yb       [auto]     optional, -yBound - +yBound of insulation plot
	--bg       []         FLAG, use transparent background of insulation plot
	--minDist  []         minimum allowed interaction distance, exclude < N distance (in BP)
	--maxDist  []         maximum allowed interaction distance, exclude > N distance (in BP)
	--ez       []         FLAG, ignore 0s in all calculations

Notes:
    This script calculates the insulation index of a given matrix to identify TAD boundaries.
    Matrix can be TXT or gzipped TXT.
    See website for matrix format details.

Contact:
    Bryan R. Lajoie
    Dekker Lab 2015
    https://github.com/blajoie/cworld-dekker
    http://my5C.umassmed.edu
```


```
Tool:		matrix2insulationRange.pl
Version:	0.01
Summary:	calculate insulation index (TADs) over range of square sizes

Usage: perl matrix2insulationRange.pl [OPTIONS] -i <inputMatrix>

Required:
	-i         []         input matrix file

Options:
	-v         []         FLAG, verbose mode
	-o         []         prefix for output file(s)
	--im       []         optional, insulation mode (how to aggregrate signal within insulation square), mean,sum,median
	--istart   []         insulation square size start, minumum bin size
	--iend     []         insulation square size end, maximum bin size
	--istep    []         insulation square size step, step size (in bp) for insulation square range
	--minDist  []         minimum allowed interaction distance, exclude < N distance (in BP)
	--maxDist  []         maximum allowed interaction distance, exclude > N distance (in BP)
	--ez       []         FLAG, ignore 0s in all calculations

Notes:
    This script calculates the insulation index over a range of insulation square sizes.
    Matrix can be TXT or gzipped TXT.
    See website for matrix format details.

Contact:
    Bryan R. Lajoie
    Dekker Lab 2015
    https://github.com/blajoie/cworld-dekker
    http://my5C.umassmed.edu
```


```
Tool:		matrix2loess.pl
Version:	0.01
Summary:	calculate the loess (expected/stdev/zScore) for a given matrix

Usage: perl matrix2loess.pl [OPTIONS] -i <inputMatrix>

Required:
	-i         []         input matrix file

Options:
	-v         []         FLAG, verbose mode
	-o         []         prefix for output file(s)
	--lof      []         optional loess object file (pre-calculated loess)
	--ec       []         FLAG, exclude CIS data
	--et       []         FLAG, exclude TRANS data
	--ez       []         FLAG, ignore zeros from all calculations
	--ca       [0.01]     lowess alpha value, fraction of datapoints to smooth over
	--caf      [1000]     cis approximate factor to speed up loess, genomic distance / -caffs
	--dif      []         FLAG, disable loess IQR (outlier) filter
	--minDist  []         minimum allowed interaction distance, exclude < N distance (in BP)
	--maxDist  []         maximum allowed interaction distance, exclude > N distance (in BP)
	--ez       []         FLAG, ignore 0s in all calculations
	--smf      []         FLAG, supress all output matrix files
	--slf      []         FLAG, supress loess file
	--d        []         FLAG, enable debug mode

Notes:
    This script calculates a lowess smoothing of the data per distance.
    Matrix can be TXT or gzipped TXT.
    See website for matrix format details.

Contact:
    Bryan R. Lajoie
    Dekker Lab 2015
    https://github.com/blajoie/cworld-dekker
    http://my5C.umassmed.edu
```


```
Tool:		matrix2pairwise.pl
Version:	0.01
Summary:	transform tsv matrix into 3 column tsv file

Usage: perl matrix2pairwise.pl [OPTIONS] -i <inputMatrix>

Required:
	-i         []         input matrix file

Options:
	-v         []         FLAG, verbose mode
	-o         []         prefix for output file(s)
	--ec       []         FLAG, exclude CIS data
	--et       []         FLAG, exclude TRANS data
	--sna      []         FLAG, skip NAs in output
	--ez       []         FLAG, ignore zeros from all calculations

Notes:
    This script calculates a lowess smoothing of the data per distance.
    Matrix can be TXT or gzipped TXT.
    See website for matrix format details.

Contact:
    Bryan R. Lajoie
    Dekker Lab 2015
    https://github.com/blajoie/cworld-dekker
    http://my5C.umassmed.edu
```


```
Tool:		matrix2scaling.pl
Version:	0.01
Summary:	transform matrix into scaling (polymer) plot

Usage: perl matrix2scaling.pl [OPTIONS] -i <inputMatrix>

Required:
	-i         []         input matrix file (can accept multiple files)

Options:
	-v         []         FLAG, verbose mode
	-o         []         prefix for output file(s)
	--lof      []         optional loess object file (pre-calculated loess)
	--n        []         optional job name for output files
	--ca       [0.01]     lowess alpha value, fraction of datapoints to smooth over
	--dif      []         FLAG, disable loess IQR (outlier) filter
	--minDist  []         minimum allowed interaction distance, exclude < N distance (in BP)
	--maxDist  []         maximum allowed interaction distance, exclude > N distance (in BP)
	--caf      [1000]     cis approximate factor to speed up loess, genomic distance / -caffs
	--ez       []         FLAG, ignore 0s in all calculations

Notes:
    This script calculates a scaling plot of a given matrix.
    Matrix can be TXT or gzipped TXT.
    See website for matrix format details.

Contact:
    Bryan R. Lajoie
    Dekker Lab 2015
    https://github.com/blajoie/cworld-dekker
    http://my5C.umassmed.edu
```


```
Tool:		matrix2stacked.pl
Version:	0.01
Summary:	transform matrix into stacked anchor matrix

Usage: perl matrix2stacked.pl [OPTIONS] -i <inputMatrix>

Required:
	-i         []         input matrix file (can accept multiple files)

Options:
	-v         []         FLAG, verbose mode
	-o         []         prefix for output file(s)
	--minDist  []         minimum allowed interaction distance, exclude < N distance (in BP)
	--maxDist  []         maximum allowed interaction distance, exclude > N distance (in BP)

Notes:
    This script transforms a matrix into stacked anchor matrix.
    Matrix can be TXT or gzipped TXT.
    See website for matrix format details.

Contact:
    Bryan R. Lajoie
    Dekker Lab 2015
    https://github.com/blajoie/cworld-dekker
    http://my5C.umassmed.edu
```


```
Tool:		matrix2symmetrical.pl
Version:	0.01
Summary:	transform rectangular matrix into square (symmetrical) matrix

Usage: perl matrix2symmetrical.pl [OPTIONS] -i <inputMatrix>

Required:
	-i         []         input matrix file

Options:
	-v         []         FLAG, verbose mode
	-o         []         prefix for output file(s)

Notes:
    This script turns any matrix into a symmetrical matrix.  If the matrix is already symmetrical, nothing is changed.
    Matrix can be TXT or gzipped TXT.
    See website for matrix format details.

Contact:
    Bryan R. Lajoie
    Dekker Lab 2015
    https://github.com/blajoie/cworld-dekker
    http://my5C.umassmed.edu
```


```
Tool:		matrix2webplot.pl
Version:	0.01
Summary:	draws 'web-plot' of matrix file

Usage: perl matrix2webplot.pl [OPTIONS] -i <inputMatrix>

Required:
	-i@        []         input matrix file [MULTIPLE]

Options:
	-v         []         FLAG, verbose mode
	-o         []         prefix for output file(s)
	--minDist  []         minimum allowed interaction distance, exclude < N distance (in BP)
	--maxDist  []         maximum allowed interaction distance, exclude > N distance (in BP)
	--bg       []         FLAG, use transparent background
	--maxdim   [16000]    maximum image dimension [hard-limited:32000]
	--lt       [0]        log transform input data into specified base
	--start    []         absolute value for color start
	--end      []         absolute vlaue for color end
	--startTile [0.025]    fraction value for color start
	--endTile  [0.975]    fraction value for color end
	--pc       [white,orange,red,darkRed] positive color scale string
	--nc       [white,cyan,blue,darkBlue] negative color scale string
	--mc       [null]     missing data color

Notes:
    This script coverts a matrix into a PNG matrix2webplot image.
    Matrix should have row/col headers in the standard my5C format.
    Matrix can be TXT or gzipped TXT.
    See website for matrix format details.

Contact:
    Bryan R. Lajoie
    Dekker Lab 2015
    https://github.com/blajoie/cworld-dekker
    http://my5C.umassmed.edu
```


```
Tool:		normalizeMatrix.pl
Version:	0.01
Summary:	normalizes matrix sum - scales to 10^6

Usage: perl normalizeMatrix.pl [OPTIONS] -i <inputMatrix>

Required:
	-i         []         input matrix file, multiple files allowed

Options:
	-v         []         FLAG, verbose mode
	-o         []         prefix for output file(s)
	--id       []         FLAG, ignore diagonal bin during normalization

Notes:
    This script normalizes the sum of a matrix.
    If the matrix is symmetrical, then only half of the matrix is summed (diagonal + x>y data)
    Each interaction is first translated into a percentage of the current sum.
    Each value is then multipled by 1,000,000 so that the new sum is 10^6
    Matrix should have row/col headers in the standard my5C format.
    Matrix can be TXT or gzipped TXT.
    See website for matrix format details.

Contact:
    Bryan R. Lajoie
    Dekker Lab 2015
    https://github.com/blajoie/cworld-dekker
    http://my5C.umassmed.edu
```


```
Tool:		primer2plates.pl
Version:	0.01
Summary:	layout primers in 96-well plate format

Usage: perl primer2plates.pl [OPTIONS] -i <inputPrimers>

Required:
	-i         []         input matrix file, multiple files allowed

Options:
	-v         []         FLAG, verbose mode
	-o         []         prefix for output file(s)

Notes:
    This script can layout primers in 96-well plate format
    Matrix should have row/col headers in the standard my5C format.
    Matrix can be TXT or gzipped TXT.
    See website for matrix format details.

Contact:
    Bryan R. Lajoie
    Dekker Lab 2015
    https://github.com/blajoie/cworld-dekker
    http://my5C.umassmed.edu
```


```
Tool:		reOrderMatrix.pl
Version:	0.01
Summary:	re-order matrix by list of headers

Usage: perl reOrderMatrix.pl [OPTIONS] -i <inputMatrix>

Required:
	-i         []         input matrix file
	--xohl     []         x-axis headers list file
	--yohl     []         y-axis headers list file

Options:
	-v         []         FLAG, verbose mode
	-o         []         prefix for output file(s)

Notes:
    This script can layout primers in 96-well plate format
    Matrix should have row/col headers in the standard my5C format.
    Matrix can be TXT or gzipped TXT.
    See website for matrix format details.

Contact:
    Bryan R. Lajoie
    Dekker Lab 2015
    https://github.com/blajoie/cworld-dekker
    http://my5C.umassmed.edu
```


```
Tool:		singletonRemoval.pl
Version:	0.01
Summary:	detect and remove singleton outliers

Usage: perl singletonRemoval.pl [OPTIONS] -i <inputMatrix>

Required:
	-i         []         input matrix file

Options:
	-v         []         FLAG, verbose mode
	-o         []         prefix for output file(s)
	--lof      []         optional loess object file (pre-calculated loess)
	--ic       []         model CIS data to detect outlier row/cols
	--it       []         model TRANS data to detect outlier row/cols
	--mof      []         optional manual outlier file, 1 header per line to be filtered.
	--ca       [0.01]     lowess alpha value, fraction of datapoints to smooth over
	--dif      []         FLAG, disable loess IQR (outlier) filter
	--minDist  []         minimum allowed interaction distance, exclude < N distance (in BP)
	--maxDist  []         maximum allowed interaction distance, exclude > N distance (in BP)
	--caf      [1000]     cis approximate factor to speed up loess, genomic distance / -caffs
	--cta      [12]       z-score threshold for cis interactions
	--tta      [12]       z-score threshold for trans interactions
	--ez       []         FLAG, ignore 0s in all calculations

Notes:
    This script can detect/remove singleton outlier interactions.
    Matrix should have row/col headers in the standard my5C format.
    Matrix can be TXT or gzipped TXT.
    See website for matrix format details.

Contact:
    Bryan R. Lajoie
    Dekker Lab 2015
    https://github.com/blajoie/cworld-dekker
    http://my5C.umassmed.edu
```


```
Tool:		subsetMatrix.pl
Version:	0.01
Summary:	subset matrix by distance, or by BED file (bin overlap)

Usage: perl subsetMatrix.pl [OPTIONS] -i <inputMatrix>

Required:
	-i         []         input matrix file

Options:
	-v         []         FLAG, verbose mode
	-o         []         prefix for output file(s)
	--minDist  []         minimum allowed interaction distance, exclude < N distance (in BP)
	--maxDist  []         maximum allowed interaction distance, exclude > N distance (in BP)
	--ec       []         FLAG, exclude CIS data
	--et       []         FLAG, exclude TRANS data
	--ebf@     []         x/y axis element bed file
	--yebf@    []         y axis element bed file
	--xebf@    []         x axis element bed file
	--z@       []         x/y axis zoom coordinate [UCSC]
	--yz@      []         y axis zoom coordinate [UCSC]
	--xz@      []         x axis zoom coordinate [UCSC]
	--en       []         elementName, descriptor for output files
	-zs        []         increase element size by N bp, (increase overlap)

Notes:
    This script aggregrates cData surrounding a list of 'elements'.
    Input Matrix can be TXT or gzipped TXT.
    AnchorListFile must be BED5+ format.
    See website for matrix format details.

Contact:
    Bryan R. Lajoie
    Dekker Lab 2015
    https://github.com/blajoie/cworld-dekker
    http://my5C.umassmed.edu
```


```
Tool:		symmetrical2seperate.pl
Version:	0.01
Summary:	transform symmetrical matrix into axis-seperated [5C headers only]

Usage: perl symmetrical2seperate.pl [OPTIONS] -i <inputMatrix>

Required:
	-i         []         input matrix file

Options:
	-v         []         FLAG, verbose mode
	-o         []         prefix for output file(s)

Notes:
    This script turns any symmetrical matrix into a axis-seperated matrix.  FORWARD on y-axis, REVERSE on x-axis.
    Matrix can be TXT or gzipped TXT.
    See website for matrix format details.

Contact:
    Bryan R. Lajoie
    Dekker Lab 2015
    https://github.com/blajoie/cworld-dekker
    http://my5C.umassmed.edu
```


```
Tool:		elementPileUp.pl
Version:	0.01
Summary:	pile up cData around specified list of 'elements'

Usage: perl elementPileUp.pl [OPTIONS] -i <inputMatrix>

Required:
	-i         []         input matrix file
	--ebf      []         element bed file

Options:
	-v         []         FLAG, verbose mode
	-o         []         prefix for output file(s)
	-zs        []         elementZoneSize, size of zone to use around element (x axis - in BP)
	--bs       [1]        bucketSpan, bucket size for discrete binning
	--ec       []         FLAG, exclude CIS data
	--et       []         FLAG, exclude TRANS data
	--ez       []         FLAG, ignore zeros from all calculations

Notes:
    This script aggregrates cData surrounding a list of 'elements'.
    Input Matrix can be TXT or gzipped TXT.
    AnchorListFile must be BED5+ format.
    See website for matrix format details.

Contact:
    Bryan R. Lajoie
    Dekker Lab 2015
    https://github.com/blajoie/cworld-dekker
    http://my5C.umassmed.edu
```
  
## Usage Examples

```
perl ../scripts/perl/addMatrixHeaders.pl -i ../sample-data/addMatrixHeaders/5C.naked.matrix.gz --xhf ../sample-data/addMatrixHeaders/5C.xHeaders.gz --yhf ../sample-data/addMatrixHeaders/5C.yHeaders.gz
perl ../scripts/perl/aggregateBED.pl -i ../sample-data/aggregrateBED/chrX.bed.gz --a mm9 --wsize 40000 --wstep 1 --wmode sum
perl ../scripts/perl/anchorPileUp.pl -i ../sample-data/anchorPileUp/NPC_chrX.matrix.gz --abf ../sample-data/anchorPileUp/Xi-escapees.bed --azs 50000000 --ez --maxDist 50000000 
perl ../scripts/perl/anchorPurge.pl -i ../sample-data/addMatrixHeaders/K5.matrix.gz --ic --ta 0.5
perl ../scripts/perl/applyCorrection.pl -i ../sample-data/addMatrixHeaders/K5.matrix.gz --ff ../sample-data/applyCorrection/K5.allPrimerFactors --ic
perl ../scripts/perl/binMatrix.pl -i ../sample-data/addMatrixHeaders/K5.matrix.gz --bsize 30000 --bstep 10 --bmode median
perl ../scripts/perl/changeMatrixHeaders.pl -i ../sample-data/addMatrixHeaders/K5.matrix.gz --hmf ../sample-data/changeMatrixHeaders/K5.headers.map 
perl ../scripts/perl/collapseMatrix.pl -i ../sample-data/collapseMatrix/NPC_chr14-chr15-chr16.matrix.gz
perl ../scripts/perl/column2matrix.pl -i ../sample-data/column2matrix/K5.pairwise.txt.gz --oxh ../sample-data/column2matrix/K5.x.headers --oyh ../sample-data/column2matrix/K5.y.headers
perl ../scripts/perl/combineMatrices.pl -i ../sample-data/addMatrixHeaders/K5.matrix.gz -i ../sample-data/combineMatrices/GM.matrix.gz --cm max -o K5__GM 
perl ../scripts/perl/compareInsulation.pl -1 ../sample-data/compareInsulation/N2.insulation.gz -2 ../sample-data/compareInsulation/SRy93.insulation.gz
perl ../scripts/perl/compareMatrices.pl -1 ../sample-data/addMatrixHeaders/K5.matrix.gz -2 ../sample-data/combineMatrices/GM.matrix.gz -o K5__GM --cm subtract
perl ../scripts/perl/correlateMatrices.pl -1 ../sample-data/addMatrixHeaders/K5.matrix.gz -2 ../sample-data/combineMatrices/GM.matrix.gz -o K5__GM
perl ../scripts/perl/coverageCorrect.pl -i ../sample-data/addMatrixHeaders/K5.matrix.gz --ic --ct 0.2
perl ../scripts/perl/digitizePicture.pl -i ../sample-data/digitizePicture/K5.png
perl ../scripts/perl/extractSubMatrices.pl -i ../sample-data/collapseMatrix/NPC_chr14-chr15-chr16.matrix.gz --ozc chr14:50000000--70000000 --ozc chr15:1--30000000 --ozc chr16:30000000--90000000 --ozc chr14:1--20000000
perl ../scripts/perl/fillMissingData.pl -i ../sample-data/fillMissingData/K5.outlierFiltered.matrix.gz
perl ../scripts/perl/generateBins.pl --r chr1:1-2000000 --bsize 30000 --bstep 10 -a hg19 chr1-FGF5
perl ../scripts/perl/heatmap.pl -i ../sample-data/binMatrix/K5__30000__10.matrix.gz -i ../sample-data/binMatrix/GM__30000__10.matrix.gz -o K5__GM__30000__10__median
perl ../scripts/perl/heatmap.pl -i ../sample-data/anchorPileUp/NPC_chrX.matrix.gz --hbf ../sample-data/anchorPileUp/Xi-escapees.bed --hc cyan
perl ../scripts/perl/insulation2tads.pl -i ../sample-data/insulation2tads/NPC.insulation.gz --b ../sample-data/insulation2tads/NPC.boundaries.gz
perl ../scripts/perl/interactionPileUp.pl -i ../sample-data/interactionPileUp/MeyersHiC-N2-DpnII.ce10.NA.H-10000-wDiag-noSS-iced.normalized.subset4000000__chrX__chrX__cis.matrix.gz --ebf ../sample-data/interactionPileUp/Top25RexPRex.bed --ezs 200000 --am sum --maxED 3
perl ../scripts/perl/matrix2anchorPlot.pl -i ../sample-data/addMatrixHeaders/K5.matrix.gz
perl ../scripts/perl/matrix2bed12.pl -i ../sample-data/matrix2anchorPlot/K5-highlight.matrix --start 1
perl ../scripts/perl/matrix2compartment.pl -i ../sample-data/collapseMatrix/NPC_chr14-chr15-chr16.matrix.gz
perl ../scripts/perl/matrix2compartment.pl -i ../sample-data/collapseMatrix/NPC_chr14-chr15-chr16.matrix.gz --ec
perl ../scripts/perl/matrix2direction.pl -i ../sample-data/collapseMatrix/NPC_chr14-chr15-chr16__chr14__chr14__cis.matrix.gz --ds 5000000
perl ../scripts/perl/matrix2distance.pl -i ../sample-data/addMatrixHeaders/K5.matrix.gz --sn 
perl ../scripts/perl/matrix2headerBed.pl -i ../sample-data/collapseMatrix/NPC_chr14-chr15-chr16.matrix.gz 
perl ../scripts/perl/matrix2info.pl -i ../sample-data/addMatrixHeaders/K5.matrix.gz
perl ../scripts/perl/matrix2insulationRange.pl -i ../sample-data/binMatrix/K5__30000__10.matrix.gz --istart 40000 --iend 4000000 --istep 40000 --im mean
perl ../scripts/perl/matrix2insulation.pl -i ../sample-data/collapseMatrix/NPC_chr14-chr15-chr16__chr14__chr14__cis.matrix.gz --is 480000 --ids 320000 --im iqrMean --nt 0 --ss 160000 --yb 1.5 --nt 0 --bmoe 0
perl ../scripts/perl/matrix2loess.pl -i ../sample-data/addMatrixHeaders/K5.matrix.gz --et
perl ../scripts/perl/matrix2pairwise.pl -i ../sample-data/addMatrixHeaders/K5.matrix.gz --et --sna --ez
perl ../scripts/perl/matrix2scaling.pl -i ../sample-data/addMatrixHeaders/K5.matrix.gz -i ../sample-data/combineMatrices/GM.matrix.gz -o K5__GM
perl ../scripts/perl/matrix2stacked.pl -i ../sample-data/collapseMatrix/NPC_chr14-chr15-chr16__chr14__chr14__cis.matrix.gz --maxDist 5000000
perl ../scripts/perl/matrix2symmetrical.pl -i ../sample-data/addMatrixHeaders/K5.matrix.gz
perl ../scripts/perl/matrix2webplot.pl -i ../sample-data/matrix2symmetrical/K5-highlight.symmetrical.matrix.gz
perl ../scripts/perl/normalizeMatrix.pl -i ../sample-data/addMatrixHeaders/K5.matrix.gz -i ../sample-data/combineMatrices/GM.matrix.gz
perl ../scripts/perl/reOrderMatrix.pl -i ../sample-data/addMatrixHeaders/K5.matrix.gz --xohl ../sample-data/reOrderMatrix/K5.rev --yohl ../sample-data/reOrderMatrix/K5.for
perl ../scripts/perl/singletonRemoval.pl -i ../sample-data/collapseMatrix/NPC_chr14-chr15-chr16.matrix.gz --it --ic --cta 2 --tta 2
perl ../scripts/perl/subsetMatrix.pl -i ../sample-data/collapseMatrix/NPC_chr14-chr15-chr16.matrix.gz -z chr14:20000000-40000000 --yz chr15:4000000-100000000 --xz chr16:20000000-100000000
perl ../scripts/perl/subsetMatrix.pl -i ../sample-data/collapseMatrix/NPC_chr14-chr15-chr16.matrix.gz --yz chr14:30000000-60000000 --yz chr15:20000000-50000000 --yz chr16:30000000-40000000 -z chr14:1-200000000
perl ../scripts/perl/symmetrical2seperate.pl -i ../sample-data/symmetrical2seperate/K5.symmetrical.matrix.gz
perl ../scripts/perl/tickPlot.pl -i ../sample-data/collapseMatrix/NPC_chr14-chr15-chr16__chr14__chr14__cis.matrix.gz --ebf ../sample-data/tickPlot/chr14.bed --bs 25
```

## Change Log

## Bugs and Feedback

For bugs, questions and discussions please use the [Github Issues](https://github.com/blajoie/cworld::dekker/issues).

## LICENSE

Licensed under the Apache License, Version 2.0 (the 'License');
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

<http://www.apache.org/licenses/LICENSE-2.0>

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an 'AS IS' BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

