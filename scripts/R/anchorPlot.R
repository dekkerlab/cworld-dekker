args <- commandArgs(TRUE)
dir<-args[1]
anchorPlotFile<-args[2]
expectedPlotFile<-args[3]
anchorName<-args[4]
anchorStart<-as.numeric(args[5])
anchorEnd<-as.numeric(args[6])
ymin<-as.numeric(args[7])
ymax<-as.numeric(args[8])

wd<-paste(dir,sep='')
setwd(wd)

expectedData<-read.table(expectedPlotFile,header=T,sep="\t",stringsAsFactors=FALSE)
anchorData<-read.table(anchorPlotFile,header=T,sep="\t",stringsAsFactors=FALSE)

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
allData.y.topIndex<-floor(allData.y.size*.995)
allData.y.bottomIndex<-ceiling(allData.y.size*0.005)
allData.y.lim.top<-allData.y[allData.y.topIndex]
allData.y.lim.bottom<-allData.y[allData.y.bottomIndex]

if(is.na(ymin)) {
	ymin<-allData.y.lim.bottom
}
if(is.na(ymax)) {
	ymax<-allData.y.lim.top
}

ymin
ymax

# now plot anchorData
pngfile<-paste(anchorPlotFile,".png",sep='')
png(pngfile,height=400,width=1200)

par(mar=c(4, 4, 4, 4))

# plot all the data - setup area
plot(allData.x,allData.y,ylim=c(ymin,ymax),xlab="Genomic Position (bp)",ylab="C counts",main=paste("Anchor Plot - ",anchorName,sep=""),type="n")

# draw anchor box
rect(anchorStart, ymin, anchorEnd, ymax, col="orange", lwd=5, border="orange")

peakData<-subset(anchorData,anchorData$highlight==1)

anchorData<-anchorData[with(anchorData, order(interactionMidpoint)), ]

# segment the anchorData data into upstream/downstream
anchorData.upstream <- subset(anchorData,(anchorData$interactionMidpoint <= anchorMidpoint))
anchorData.downstream <- subset(anchorData,(anchorData$interactionMidpoint >= anchorMidpoint))	

# cData score
lines(anchorData.upstream$interactionMidpoint,anchorData.upstream$cScore,col="red",lwd=2)
lines(anchorData.downstream$interactionMidpoint,anchorData.downstream$cScore,col="red",lwd=2)

#draw final iteration points
points(anchorData$interactionMidpoint,anchorData$cScore,col="black",pch=21,bg="red")		

#draw peak points
points(peakData$interactionMidpoint,peakData$cScore,col="black",bg=rgb(0,1,0,0.7),pch=21,cex=3)

xRange<-(region.bound.right-region.bound.left)
yRange<-(ymax-ymin)

if((nrow(peakData) > 0) & (nrow(peakData) < 20)) {
	for (i in 1:nrow(peakData) ) {
		tmpX<-peakData$interactionMidpoint[i]
		tmpY<-peakData$cScore[i]
		peakName<-peakData$xHeader[i]
		tmpPeakNameList<-unlist(strsplit(peakName, "\\|"))
		peakName<-tmpPeakNameList[1]
		text((tmpX+(xRange/100)),(tmpY+(yRange/100)),labels=peakName,adj=0,xpd=TRUE,cex=0.80,srt=35)
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

dev.off()