args <- commandArgs(TRUE)
dir<-args[1]
anchorPlotFile<-args[2]
expectedPlotFile<-args[3]
anchorName<-args[4]
anchorStart<-as.numeric(args[5])
anchorEnd<-as.numeric(args[6])
primerPerformanceFactor<-as.numeric(args[7])

wd<-paste(dir,sep='')
setwd(wd)

expectedData<-read.table(expectedPlotFile,header=T,sep="\t")
anchorData<-read.table(anchorPlotFile,header=T,sep="\t")

# adjust the anchor coordinates to be at least 0.5% width of span
anchorMidpoint<-round((anchorStart+anchorEnd)/2)
allData.range<-(max(anchorData$interactionMidpoint)-min(anchorData$interactionMidpoint))
allData.minspan<-round(allData.range*0.005)
anchorStart <- anchorMidpoint - allData.minspan
anchorEnd <- anchorMidpoint + allData.minspan

# region bounds
region.bound.left<-min(anchorData$interactionMidpoint)
region.bound.right<-max(anchorData$interactionMidpoint)

# subset the expected data to only include regions bounds
expectedData <- subset(expectedData,(expectedData$interactionMidpoint >= region.bound.left & expectedData$interactionMidpoint <= region.bound.right))

# sort the expected data
expectedData<-expectedData[with(expectedData, order(interactionMidpoint)), ]	

# segment the expected data into upstream/downstream
expectedData.upstream <- subset(expectedData,(expectedData$interactionMidpoint <= anchorStart))
expectedData.downstream <- subset(expectedData,(expectedData$interactionMidpoint >= anchorEnd))

# calculate trimmed bounds
allData.y<-c(anchorData$cScore,expectedData$loessExpected,expectedData$loessExpected+expectedData$loessStdev,expectedData$loessExpected-expectedData$loessStdev)
allData.x<-c(anchorData$interactionMidpoint,expectedData$interactionMidpoint,expectedData$interactionMidpoint,expectedData$interactionMidpoint)
# sort and calculate bounds
allData.y<-sort(allData.y)
allData.y.size<-length(allData.y)
allData.y.topIndex<-floor(allData.y.size*0.995)-1
allData.y.bottomIndex<-ceiling(allData.y.size*0.005)+1
allData.y.lim.top<-allData.y[allData.y.topIndex]
allData.y.lim.bottom<-allData.y[allData.y.bottomIndex]


distinctGroups<-unique(anchorData$iterationNumber)
distinctColors<-rainbow(length(distinctGroups))

# now plot anchorData
pngfile<-paste(anchorPlotFile,".png",sep='')
png(pngfile,height=300,width=900)

if(length(distinctGroups) > 1) {
	par(mar=c(5.1, 4.1, 4.1, 8.1))
}

# plot all the data - setup area
plot(allData.x,allData.y,ylim=c(allData.y.lim.bottom,allData.y.lim.top),xlab="Genomic Position (bp)",ylab="cScore",main=paste(anchorName,"\n","Anchor Plot","\n","primer factor = ",primerPerformanceFactor,sep=""),type="n")

# draw anchor box
rect(anchorStart, allData.y.lim.bottom, anchorEnd, allData.y.lim.top, col="orange", lwd=5, border="orange")


for(group in distinctGroups) {
	
	if(group == max(distinctGroups)) {
		lineWidth=2
	} else {
		lineWidth=1
	}
	
	tmp <- subset(anchorData,anchorData$iterationNumber==group)
	tmp<-tmp[with(tmp, order(interactionMidpoint)), ]

	# segment the tmp data into upstream/downstream
	tmp.upstream <- subset(tmp,(tmp$interactionMidpoint <= anchorMidpoint))
	tmp.downstream <- subset(tmp,(tmp$interactionMidpoint >= anchorMidpoint))	
	
	# cData score
	lines(tmp.upstream$interactionMidpoint,tmp.upstream$cScore,col=distinctColors[group+1],lwd=lineWidth)
	lines(tmp.downstream$interactionMidpoint,tmp.downstream$cScore,col=distinctColors[group+1],lwd=lineWidth)
	
	#draw raw data points
	if(group == min(distinctGroups)) {	
		points(tmp$interactionMidpoint,tmp$cScore,col=distinctColors[group+1],cex=0.75)		
	}
	
	#draw final iteration points
	if(group == max(distinctGroups)) {	
		points(tmp$interactionMidpoint,tmp$cScore,col="black",pch=21,bg=distinctColors[group+1])		
	}
}

# plot the full expected data & +/- loess stdev lines
#upstream
lines(expectedData.upstream$interactionMidpoint,expectedData.upstream$loessExpected,col="black",lwd=2)
lines(expectedData.upstream$interactionMidpoint,expectedData.upstream$loessExpected+expectedData.upstream$loessStdev,col="black",lwd=1,lty=2)
lines(expectedData.upstream$interactionMidpoint,expectedData.upstream$loessExpected-expectedData.upstream$loessStdev,col="black",lwd=1,lty=2)
#downstream
lines(expectedData.downstream$interactionMidpoint,expectedData.downstream$loessExpected,col="black",lwd=2)
lines(expectedData.downstream$interactionMidpoint,expectedData.downstream$loessExpected+expectedData.downstream$loessStdev,col="black",lwd=1,lty=2)
lines(expectedData.downstream$interactionMidpoint,expectedData.downstream$loessExpected-expectedData.downstream$loessStdev,col="black",lwd=1,lty=2)

if(length(distinctGroups) > 1) {
	# allow legend to be plotted outside of plot area
	par(xpd=TRUE) 
	legend("topright", inset=c(-0.105,0), legend=paste("iteration # ",unique(anchorData$iterationNumber),sep=""), text.col=distinctColors,cex=0.75)
}

dev.off()