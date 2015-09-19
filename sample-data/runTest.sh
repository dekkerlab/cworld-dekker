perl -I ../lib ../scripts/perl/addMatrixHeaders.pl -i ../sample-data/addMatrixHeaders/5C.naked.matrix.gz --xhf ../sample-data/addMatrixHeaders/5C.xHeaders.gz --yhf ../sample-data/addMatrixHeaders/5C.yHeaders.gz
perl -I ../lib ../scripts/perl/aggregateBED.pl -i ../sample-data/aggregrateBED/chrX.bed.gz --a mm9 --wsize 40000 --wstep 1 --wmode sum
perl -I ../lib ../scripts/perl/elementPileUp.pl -i ../sample-data/elementPileUp/NPC_chrX.matrix.gz --ebf ../sample-data/elementPileUp/Xi-escapees.bed --ezs 50000000 --ez --maxDist 50000000 
perl -I ../lib ../scripts/perl/anchorPurge.pl -i ../sample-data/addMatrixHeaders/K5.matrix.gz --ic --ta 0.5
perl -I ../lib ../scripts/perl/applyCorrection.pl -i ../sample-data/addMatrixHeaders/K5.matrix.gz --ff ../sample-data/applyCorrection/K5.allPrimerFactors --ic
perl -I ../lib ../scripts/perl/binMatrix.pl -i ../sample-data/addMatrixHeaders/K5.matrix.gz --bsize 30000 --bstep 10 --bmode median
perl -I ../lib ../scripts/perl/changeMatrixHeaders.pl -i ../sample-data/addMatrixHeaders/K5.matrix.gz --hmf ../sample-data/changeMatrixHeaders/K5.headers.map 
perl -I ../lib ../scripts/perl/collapseMatrix.pl -i ../sample-data/collapseMatrix/NPC_chr14-chr15-chr16.matrix.gz
perl -I ../lib ../scripts/perl/column2matrix.pl -i ../sample-data/column2matrix/K5.pairwise.txt.gz --oxh ../sample-data/column2matrix/K5.x.headers --oyh ../sample-data/column2matrix/K5.y.headers
perl -I ../lib ../scripts/perl/combineMatrices.pl -i ../sample-data/addMatrixHeaders/K5.matrix.gz -i ../sample-data/combineMatrices/GM.matrix.gz --cm max -o K5__GM 
perl -I ../lib ../scripts/perl/compareInsulation.pl -1 ../sample-data/compareInsulation/N2.insulation.gz -2 ../sample-data/compareInsulation/SRy93.insulation.gz
perl -I ../lib ../scripts/perl/compareMatrices.pl -1 ../sample-data/addMatrixHeaders/K5.matrix.gz -2 ../sample-data/combineMatrices/GM.matrix.gz -o K5__GM --cm subtract
perl -I ../lib ../scripts/perl/correlateMatrices.pl -1 ../sample-data/addMatrixHeaders/K5.matrix.gz -2 ../sample-data/combineMatrices/GM.matrix.gz -o K5__GM
perl -I ../lib ../scripts/perl/coverageCorrect.pl -i ../sample-data/addMatrixHeaders/K5.matrix.gz --ic --ct 0.2
perl -I ../lib ../scripts/perl/digitizePicture.pl -i ../sample-data/digitizePicture/K5.png
perl -I ../lib ../scripts/perl/extractSubMatrices.pl -i ../sample-data/collapseMatrix/NPC_chr14-chr15-chr16.matrix.gz --ozc chr14:50000000--70000000 --ozc chr15:1--30000000 --ozc chr16:30000000--90000000 --ozc chr14:1--20000000
perl -I ../lib ../scripts/perl/fillMissingData.pl -i ../sample-data/fillMissingData/K5.outlierFiltered.matrix.gz
perl -I ../lib ../scripts/perl/generateBins.pl --r chr1:1-2000000 --bsize 30000 --bstep 10 -a hg19 chr1-FGF5
perl -I ../lib ../scripts/perl/heatmap.pl -i ../sample-data/binMatrix/K5__30000__10.matrix.gz -i ../sample-data/binMatrix/GM__30000__10.matrix.gz -o K5__GM__30000__10__median
perl -I ../lib ../scripts/perl/heatmap.pl -i ../sample-data/elementPileUp/NPC_chrX.matrix.gz --hbf ../sample-data/elementPileUp/Xi-escapees.bed --hc cyan
perl -I ../lib ../scripts/perl/insulation2tads.pl -i ../sample-data/insulation2tads/NPC.insulation.gz --b ../sample-data/insulation2tads/NPC.boundaries.gz
perl -I ../lib ../scripts/perl/interactionPileUp.pl -i ../sample-data/interactionPileUp/MeyersHiC-N2-DpnII.ce10.NA.H-10000-wDiag-noSS-iced.normalized.subset4000000__chrX__chrX__cis.matrix.gz --ebf ../sample-data/interactionPileUp/Top25RexPRex.bed --ezs 200000 --am sum --maxED 3
perl -I ../lib ../scripts/perl/matrix2anchorPlot.pl -i ../sample-data/addMatrixHeaders/K5.matrix.gz
perl -I ../lib/ ../scripts/perl/matrix2bed12.pl -i ../sample-data/matrix2anchorPlot/K5-highlight.matrix --start 1
perl -I ../lib/ ../scripts/perl/matrix2compartment.pl -i ../sample-data/collapseMatrix/NPC_chr14-chr15-chr16.matrix.gz
perl -I ../lib/ ../scripts/perl/matrix2compartment.pl -i ../sample-data/collapseMatrix/NPC_chr14-chr15-chr16.matrix.gz --ec
perl -I ../lib/ ../scripts/perl/matrix2direction.pl -i ../sample-data/collapseMatrix/NPC_chr14-chr15-chr16__chr14__chr14__cis.matrix.gz --ds 5000000
perl -I ../lib/ ../scripts/perl/matrix2distance.pl -i ../sample-data/addMatrixHeaders/K5.matrix.gz --sn 
perl -I ../lib/ ../scripts/perl/matrix2headerBed.pl -i ../sample-data/collapseMatrix/NPC_chr14-chr15-chr16.matrix.gz 
perl -I ../lib/ ../scripts/perl/matrix2info.pl -i ../sample-data/addMatrixHeaders/K5.matrix.gz
perl -I ../lib ../scripts/perl/matrix2insulationRange.pl -i ../sample-data/binMatrix/K5__30000__10.matrix.gz --istart 40000 --iend 4000000 --istep 40000 --im mean
perl -I ../lib/ ../scripts/perl/matrix2insulation.pl -i ../sample-data/collapseMatrix/NPC_chr14-chr15-chr16__chr14__chr14__cis.matrix.gz --is 480000 --ids 320000 --im iqrMean --nt 0 --ss 160000 --yb 1.5 --nt 0 --bmoe 0
perl -I ../lib ../scripts/perl/matrix2loess.pl -i ../sample-data/addMatrixHeaders/K5.matrix.gz --et
perl -I ../lib ../scripts/perl/matrix2pairwise.pl -i ../sample-data/addMatrixHeaders/K5.matrix.gz --et --sna --ez
perl -I ../lib ../scripts/perl/matrix2scaling.pl -i ../sample-data/addMatrixHeaders/K5.matrix.gz -i ../sample-data/combineMatrices/GM.matrix.gz -o K5__GM
perl -I ../lib ../scripts/perl/matrix2stacked.pl -i ../sample-data/collapseMatrix/NPC_chr14-chr15-chr16__chr14__chr14__cis.matrix.gz --maxDist 5000000
perl -I ../lib ../scripts/perl/matrix2symmetrical.pl -i ../sample-data/addMatrixHeaders/K5.matrix.gz
perl -I ../lib ../scripts/perl/matrix2webplot.pl -i ../sample-data/matrix2symmetrical/K5-highlight.symmetrical.matrix.gz
perl -I ../lib ../scripts/perl/normalizeMatrix.pl -i ../sample-data/addMatrixHeaders/K5.matrix.gz -i ../sample-data/combineMatrices/GM.matrix.gz
perl -I ../lib ../scripts/perl/reOrderMatrix.pl -i ../sample-data/addMatrixHeaders/K5.matrix.gz --xohl ../sample-data/reOrderMatrix/K5.rev --yohl ../sample-data/reOrderMatrix/K5.for
perl -I ../lib ../scripts/perl/singletonRemoval.pl -i ../sample-data/collapseMatrix/NPC_chr14-chr15-chr16.matrix.gz --it --ic --cta 2 --tta 2
perl -I ../lib ../scripts/perl/subsetMatrix.pl -i ../sample-data/collapseMatrix/NPC_chr14-chr15-chr16.matrix.gz -z chr14:20000000-40000000 --yz chr15:4000000-100000000 --xz chr16:20000000-100000000
perl -I ../lib ../scripts/perl/subsetMatrix.pl -i ../sample-data/collapseMatrix/NPC_chr14-chr15-chr16.matrix.gz --yz chr14:30000000-60000000 --yz chr15:20000000-50000000 --yz chr16:30000000-40000000 -z chr14:1-200000000
perl -I ../lib ../scripts/perl/symmetrical2seperate.pl -i ../sample-data/symmetrical2seperate/K5.symmetrical.matrix.gz
perl -I ../lib ../scripts/perl/tickPlot.pl -i ../sample-data/collapseMatrix/NPC_chr14-chr15-chr16__chr14__chr14__cis.matrix.gz --ebf ../sample-data/tickPlot/chr14.bed --bs 25