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
imagemagick - http://www.imagemagick.org/script/index.php

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

Grab the latest stable release
```
https://github.com/blajoie/cworld-dekker/releases
```

or clone the git project
```
[ssh] - git clone git@github.com:blajoie/cworld-dekker.git
[https] - git clone https://github.com/blajoie/cworld-dekker.git
```

To install the module:
```
perl Build.PL
./Build
./Build install
```
OR 
```
perl Build.PL
./Build
./Build install --install_base /your/custom/dir
(ensure /your/custom/dir is added to your PERL5LIB path)

e.g.
./Build install --install_base ~/perl5
# then in .bashrc
export PERL5LIB=${PERL5LIB}:/home/<yourusername>/perl5/lib/perl5


```

After installing the module, you should be free to run anything in scripts/
e.g.
```
$ perl scripts/heatmap.pl
```
## Install Troubleshooting

Trouble with libgd?
```
libgd 2.0.33 or higher required for copyRotated support
```

You need to install libgd 2.0.33 or higher.
Try to compile from source.

Download and install libgd using this link - http://www.boutell.com/gd/http/gd-2.0.33.tar.gz
```
$ cd gd-2.0.33
$ ./configure
$ make
$ sudo make install
```
More here - http://www.physics.buffalo.edu/phy410-505/tools/install/

* You will then need to re-build GD!

Uninstall GD
e.g.
```
$ cpanp u GD
```

re-install GD
```
$ sudo cpan
install GD
```

This looks to be a known bug/issue.
See below.
http://lists.freebsd.org/pipermail/freebsd-perl/2014-February/009206.html
https://github.com/lstein/Perl-GD/issues/14
http://www.ogris.de/oldnews.html

## Full Documentation

WIP

## Communication

- [Bryan Lajoie](https://github.com/blajoie)
- Twitter: [@my5C](https://twitter.com/my5C)

## What does it do?

cworld::dekker is a large collection of perl / python / R scripts for manipulating 3C/4C/5C/Hi-C data
    
## Usage Examples

```
perl scripts/perl/addMatrixHeaders.pl -i sample-data/addMatrixHeaders/5C.naked.matrix.gz --xhf ../sample-data/addMatrixHeaders/5C.xHeaders.gz --yhf ../sample-data/addMatrixHeaders/5C.yHeaders.gz
perl scripts/perl/aggregateBED.pl -i sample-data/aggregrateBED/chrX.bed.gz --a mm9 --wsize 40000 --wstep 1 --wmode sum
perl scripts/perl/elementPileUp.pl -i sample-data/elementPileUp/NPC_chrX.matrix.gz --ebf ../sample-data/elementPileUp/Xi-escapees.bed --ezs 50000000 --ez --maxDist 50000000 
perl scripts/perl/anchorPurge.pl -i sample-data/addMatrixHeaders/K5.matrix.gz --ic --ta 0.5
perl scripts/perl/applyCorrection.pl -i sample-data/addMatrixHeaders/K5.matrix.gz --ff ../sample-data/applyCorrection/K5.allPrimerFactors --ic
perl scripts/perl/binMatrix.pl -i sample-data/addMatrixHeaders/K5.matrix.gz --bsize 30000 --bstep 10 --bmode median
perl scripts/perl/changeMatrixHeaders.pl -i sample-data/addMatrixHeaders/K5.matrix.gz --hmf ../sample-data/changeMatrixHeaders/K5.headers.map 
perl scripts/perl/collapseMatrix.pl -i sample-data/collapseMatrix/NPC_chr14-chr15-chr16.matrix.gz
perl scripts/perl/column2matrix.pl -i sample-data/column2matrix/K5.pairwise.txt.gz --oxh ../sample-data/column2matrix/K5.x.headers --oyh ../sample-data/column2matrix/K5.y.headers
perl scripts/perl/combineMatrices.pl -i sample-data/addMatrixHeaders/K5.matrix.gz -i sample-data/combineMatrices/GM.matrix.gz --cm max -o K5__GM 
perl scripts/perl/compareInsulation.pl -1 ../sample-data/compareInsulation/N2.insulation.gz -2 ../sample-data/compareInsulation/SRy93.insulation.gz
perl scripts/perl/compareMatrices.pl -1 ../sample-data/addMatrixHeaders/K5.matrix.gz -2 ../sample-data/combineMatrices/GM.matrix.gz -o K5__GM --cm subtract
perl scripts/perl/correlateMatrices.pl -1 ../sample-data/addMatrixHeaders/K5.matrix.gz -2 ../sample-data/combineMatrices/GM.matrix.gz -o K5__GM
perl scripts/perl/coverageCorrect.pl -i sample-data/addMatrixHeaders/K5.matrix.gz --cm cis --ct 0.2
perl scripts/perl/digitizePicture.pl -i sample-data/digitizePicture/K5.png
perl scripts/perl/extractSubMatrices.pl -i sample-data/collapseMatrix/NPC_chr14-chr15-chr16.matrix.gz --ozc chr14:50000000--70000000 --ozc chr15:1--30000000 --ozc chr16:30000000--90000000 --ozc chr14:1--20000000
perl scripts/perl/fillMissingData.pl -i sample-data/fillMissingData/K5.outlierFiltered.matrix.gz
perl scripts/perl/generateBins.pl --r chr1:1-2000000 --bsize 30000 --bstep 10 -a hg19 chr1-FGF5
perl scripts/perl/heatmap.pl -i sample-data/binMatrix/K5__30000__10.matrix.gz -i sample-data/binMatrix/GM__30000__10.matrix.gz -o K5__GM__30000__10__median
perl scripts/perl/heatmap.pl -i sample-data/elementPileUp/NPC_chrX.matrix.gz --ebf ../sample-data/elementPileUp/Xi-escapees.bed --hc cyan
perl scripts/perl/insulation2tads.pl -i sample-data/insulation2tads/NPC.insulation.gz --b ../sample-data/insulation2tads/NPC.boundaries.gz
perl scripts/perl/interactionPileUp.pl -i sample-data/interactionPileUp/MeyersHiC-N2-DpnII.ce10.NA.H-10000-wDiag-noSS-iced.normalized.subset4000000__chrX__chrX__cis.matrix.gz --ebf ../sample-data/interactionPileUp/Top25RexPRex.bed --ezs 200000 --am sum --maxED 3
perl scripts/perl/matrix2anchorPlot.pl -i sample-data/addMatrixHeaders/K5.matrix.gz
perl scripts/perl/matrix2bed12.pl -i sample-data/matrix2anchorPlot/K5-highlight.matrix --start 1
perl scripts/perl/matrix2compartment.pl -i sample-data/collapseMatrix/NPC_chr14-chr15-chr16.matrix.gz
perl scripts/perl/matrix2compartment.pl -i sample-data/collapseMatrix/NPC_chr14-chr15-chr16.matrix.gz --ec
perl scripts/perl/matrix2direction.pl -i sample-data/collapseMatrix/NPC_chr14-chr15-chr16__chr14__chr14__cis.matrix.gz --ds 5000000
perl scripts/perl/matrix2distance.pl -i sample-data/addMatrixHeaders/K5.matrix.gz --sn 
perl scripts/perl/matrix2headerBed.pl -i sample-data/collapseMatrix/NPC_chr14-chr15-chr16.matrix.gz 
perl scripts/perl/matrix2info.pl -i sample-data/addMatrixHeaders/K5.matrix.gz
perl scripts/perl/matrix2insulationRange.pl -i sample-data/binMatrix/K5__30000__10.matrix.gz --istart 40000 --iend 4000000 --istep 40000 --im mean
perl scripts/perl/matrix2insulation.pl -i sample-data/collapseMatrix/NPC_chr14-chr15-chr16__chr14__chr14__cis.matrix.gz --is 480000 --ids 320000 --im iqrMean --nt 0 --ss 160000 --yb 1.5 --nt 0 --bmoe 0
perl scripts/perl/matrix2loess.pl -i sample-data/addMatrixHeaders/K5.matrix.gz --et
perl scripts/perl/matrix2pairwise.pl -i sample-data/addMatrixHeaders/K5.matrix.gz --et --sna --ez
perl scripts/perl/matrix2scaling.pl -i sample-data/addMatrixHeaders/K5.matrix.gz -i sample-data/combineMatrices/GM.matrix.gz -o K5__GM
perl scripts/perl/matrix2stacked.pl -i sample-data/collapseMatrix/NPC_chr14-chr15-chr16__chr14__chr14__cis.matrix.gz --maxDist 5000000
perl scripts/perl/matrix2symmetrical.pl -i sample-data/addMatrixHeaders/K5.matrix.gz
perl scripts/perl/matrix2webplot.pl -i sample-data/matrix2symmetrical/K5-highlight.symmetrical.matrix.gz
perl scripts/perl/normalizeMatrix.pl -i sample-data/addMatrixHeaders/K5.matrix.gz -i sample-data/combineMatrices/GM.matrix.gz
perl scripts/perl/reOrderMatrix.pl -i sample-data/addMatrixHeaders/K5.matrix.gz --xohl ../sample-data/reOrderMatrix/K5.rev --yohl ../sample-data/reOrderMatrix/K5.for
perl scripts/perl/singletonRemoval.pl -i sample-data/collapseMatrix/NPC_chr14-chr15-chr16.matrix.gz --it --ic --cta 2 --tta 2
perl scripts/perl/subsetMatrix.pl -i sample-data/collapseMatrix/NPC_chr14-chr15-chr16.matrix.gz -z chr14:20000000-40000000 --yz chr15:4000000-100000000 --xz chr16:20000000-100000000
perl scripts/perl/subsetMatrix.pl -i sample-data/collapseMatrix/NPC_chr14-chr15-chr16.matrix.gz --yz chr14:30000000-60000000 --yz chr15:20000000-50000000 --yz chr16:30000000-40000000 -z chr14:1-200000000
perl scripts/perl/symmetrical2seperate.pl -i sample-data/symmetrical2seperate/K5.symmetrical.matrix.gz
perl scripts/perl/tickPlot.pl -i sample-data/collapseMatrix/NPC_chr14-chr15-chr16__chr14__chr14__cis.matrix.gz --ebf ../sample-data/tickPlot/chr14.bed --bs 25
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

