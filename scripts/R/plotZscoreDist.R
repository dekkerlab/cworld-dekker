args <- commandArgs(TRUE)
dir<-args[1]
inputFiles<-args[2]
jobName<-args[3]

setwd(dir)

fileArray<-unlist(strsplit(inputFiles,","))
nFiles<-length(fileArray)

# now plot scaling data
pngfile<-paste(jobName,".png",sep='')
png(pngfile,height=900,width=900)

par(mfrow=c(3,2))

allX<-list()
allY<-list()
allData<-data.frame()
tmpSums<-list()
for (i in 1:nFiles) {
	tmpFile<-fileArray[i]

	tmpLoessData<-read.table(tmpFile,header=T,sep="\t",stringsAsFactors=FALSE)
	allData<-rbind(allData,tmpLoessData)
	
	tmpLoessData.subset<-subset(tmpLoessData,(tmpLoessData$zScore >= 0))
	tmpLoessData.subset.sum<-sum(tmpLoessData.subset$zScore)	
	tmpSums<-c(tmpSums,tmpLoessData.subset.sum)
}

cumlativeMax<-max(unlist(tmpSums))

# cut out bottom/top 1% of data (left edge)
allData.y<-sort(unique(allData$zScore))
allData.y.size<-length(allData.y)
allData.y.topIndex<-floor(allData.y.size*0.998)
allData.y.bottomIndex<-ceiling(allData.y.size*0.002)
allData.y.lim.top<-allData.y[allData.y.topIndex]
allData.y.lim.bottom<-allData.y[allData.y.bottomIndex]

allData.subset<-subset(allData,(allData$zScore >= allData.y.lim.bottom & allData$zScore <= allData.y.lim.top))

allX<-c(allData.subset$realInteractionDistance)
allY<-c(allData.subset$zScore)

pointColors<-rainbow(nFiles,alpha=0.05)
lineColors<-rainbow(nFiles,alpha=0.5)

# plot all the data - setup area
plot(allX,allY,xlab="Interaction Distance",ylab="z-score",ylim=c(allData.y.lim.bottom,allData.y.lim.top),main=paste("CData Scaling Plot - ",jobName,sep=""),type="n")
for (i in 1:nFiles) {
	tmpFile<-fileArray[i]

	tmpLoessData<-read.table(tmpFile,header=T,sep="\t",stringsAsFactors=FALSE)
	tmpLoessData<-tmpLoessData[with(tmpLoessData, order(realInteractionDistance)), ]	
	
	tmpLoessData.subset<-subset(tmpLoessData,(tmpLoessData$zScore >= allData.y.lim.bottom & tmpLoessData$zScore <= allData.y.lim.top))
	
	tmpLoessData.subset.lowess<-lowess(tmpLoessData.subset$realInteractionDistance, tmpLoessData.subset$zScore, f = 0.25, iter = 3)
	
	points(tmpLoessData.subset$realInteractionDistance,tmpLoessData.subset$zScore,col=pointColors[i],cex=0.2)	
	lines(tmpLoessData.subset.lowess,col=lineColors[i],lwd=8)
	
}

legend("top", inset=c(-0.105,0), legend=paste(fileArray,sep=""), text.col=lineColors,cex=0.5)


# plot all the data - setup area
allData.y.lim.bottom<-0
plot(allX,allY,xlab="Interaction Distance",ylab="z-score",ylim=c(allData.y.lim.bottom,allData.y.lim.top),main=paste("CData Scaling Plot - ",jobName,sep=""),type="n")
for (i in 1:nFiles) {
	tmpFile<-fileArray[i]

	tmpLoessData<-read.table(tmpFile,header=T,sep="\t",stringsAsFactors=FALSE)
	tmpLoessData<-tmpLoessData[with(tmpLoessData, order(realInteractionDistance)), ]	
	
	tmpLoessData.subset<-subset(tmpLoessData,(tmpLoessData$zScore >= allData.y.lim.bottom & tmpLoessData$zScore <= allData.y.lim.top))
	
	tmpLoessData.subset.lowess<-lowess(tmpLoessData.subset$realInteractionDistance, tmpLoessData.subset$zScore, f = 0.25, iter = 3)
	
	points(tmpLoessData.subset$realInteractionDistance,tmpLoessData.subset$zScore,col=pointColors[i],cex=0.2)	
	lines(tmpLoessData.subset.lowess,col=lineColors[i],lwd=4)
	
}

legend("top", inset=c(-0.105,0), legend=paste(fileArray,sep=""), text.col=lineColors,cex=0.5)


#
# relative cumlative sums
#

# plot all the data - setup area
allData.y.lim.bottom<-0
plot(allX,allY,xlab="Interaction Distance",ylab="cumlative z-scores",ylim=c(0,100),main=paste("CData Scaling Plot (all (+) z-scores)- ",jobName,sep=""),type="n")
for (i in 1:nFiles) {
	tmpFile<-fileArray[i]

	tmpLoessData<-read.table(tmpFile,header=T,sep="\t",stringsAsFactors=FALSE)
	tmpLoessData<-tmpLoessData[with(tmpLoessData, order(realInteractionDistance)), ]	
	
	tmpLoessData.subset<-subset(tmpLoessData,(tmpLoessData$zScore >= allData.y.lim.bottom & tmpLoessData$zScore <= allData.y.lim.top))
	
	tmp.cumsum<-cumsum(((cumsum(tmpLoessData.subset$zScore)/sum(cumsum(tmpLoessData.subset$zScore)))*100))
	lines(tmpLoessData.subset$realInteractionDistance,tmp.cumsum,col=lineColors[i],lwd=4)
	
}
legend("top", inset=c(-0.105,0), legend=paste(fileArray,sep=""), text.col=lineColors,cex=0.5)

# plot all the data - setup area
plot(allX,allY,xlab="Interaction Distance",ylab="cumlative z-scores",ylim=c(0,100),main=paste("CData Scaling Plot (>= 2 z-scores)- ",jobName,sep=""),type="n")
allData.y.lim.bottom<-2
for (i in 1:nFiles) {
	tmpFile<-fileArray[i]

	tmpLoessData<-read.table(tmpFile,header=T,sep="\t",stringsAsFactors=FALSE)
	tmpLoessData<-tmpLoessData[with(tmpLoessData, order(realInteractionDistance)), ]	
	
	tmpLoessData.subset<-subset(tmpLoessData,(tmpLoessData$zScore >= allData.y.lim.bottom & tmpLoessData$zScore <= allData.y.lim.top))
	
	tmp.cumsum<-cumsum(((cumsum(tmpLoessData.subset$zScore)/sum(cumsum(tmpLoessData.subset$zScore)))*100))
	lines(tmpLoessData.subset$realInteractionDistance,tmp.cumsum,col=lineColors[i],lwd=4)
	
}
legend("top", inset=c(-0.105,0), legend=paste(fileArray,sep=""), text.col=lineColors,cex=0.5)

#
# absolute cumlative sums
#

# plot all the data - setup area
allData.y.lim.bottom<-0
plot(allX,allY,xlab="Interaction Distance",ylab="cumlative z-scores",ylim=c(0,cumlativeMax),main=paste("CData Scaling Plot (all (+) z-scores)- ",jobName,sep=""),type="n")
for (i in 1:nFiles) {
	tmpFile<-fileArray[i]

	tmpLoessData<-read.table(tmpFile,header=T,sep="\t",stringsAsFactors=FALSE)
	tmpLoessData<-tmpLoessData[with(tmpLoessData, order(realInteractionDistance)), ]	
	
	tmpLoessData.subset<-subset(tmpLoessData,(tmpLoessData$zScore >= allData.y.lim.bottom & tmpLoessData$zScore <= allData.y.lim.top))
	
	tmp.cumsum<-cumsum(tmpLoessData.subset$zScore)
	lines(tmpLoessData.subset$realInteractionDistance,tmp.cumsum,col=lineColors[i],lwd=4)
	
}
legend("top", inset=c(-0.105,0), legend=paste(fileArray,sep=""), text.col=lineColors,cex=0.5)

# plot all the data - setup area
plot(allX,allY,xlab="Interaction Distance",ylab="cumlative z-scores",ylim=c(0,cumlativeMax),main=paste("CData Scaling Plot (>= 2 z-scores)- ",jobName,sep=""),type="n")
allData.y.lim.bottom<-2
for (i in 1:nFiles) {
	tmpFile<-fileArray[i]

	tmpLoessData<-read.table(tmpFile,header=T,sep="\t",stringsAsFactors=FALSE)
	tmpLoessData<-tmpLoessData[with(tmpLoessData, order(realInteractionDistance)), ]	
	
	tmpLoessData.subset<-subset(tmpLoessData,(tmpLoessData$zScore >= allData.y.lim.bottom & tmpLoessData$zScore <= allData.y.lim.top))
	
	tmp.cumsum<-cumsum(tmpLoessData.subset$zScore)
	lines(tmpLoessData.subset$realInteractionDistance,tmp.cumsum,col=lineColors[i],lwd=4)
	
}
legend("top", inset=c(-0.105,0), legend=paste(fileArray,sep=""), text.col=lineColors,cex=0.5)

dev.off()