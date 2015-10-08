options(bitmapType='cairo')

args <- commandArgs(TRUE)
dir<-args[1]
inputFile<-args[2]

inputFileArray<-unlist(strsplit(inputFile, "\\/"))
inputFileName<-inputFileArray[length(inputFileArray)]

wd<-paste(dir,sep='')
setwd(wd)

data<-read.table(inputFile,header=T,sep="\t")
data<-subset(data,data$loessExpectedValue!="NA")
data<-subset(data,data$loessExpectedStdev!="NA")

allData.y<-c(data$observedData,data$loessExpectedValue,data$loessExpectedValue+data$loessExpectedStdev,data$loessExpectedValue-data$loessExpectedStdev)

allData.y<-sort(allData.y)
allData.y.size<-length(allData.y)

inputFileName<-inputFile
inputFileName<-gsub(".gz", "", inputFileName)

pngfile<-paste(inputFileName,".png",sep='')
png(pngfile,height=800,width=800)

par(mfrow=c(2,1))

# first plot - all data
allData.y.topIndex<-floor(allData.y.size*0.995)-1
allData.y.bottomIndex<-ceiling(allData.y.size*0.005)+1
allData.y.lim.top<-allData.y[allData.y.topIndex]
allData.y.lim.bottom<-allData.y[allData.y.bottomIndex]

alpha<-(100/nrow(data))
if(alpha < 0.01) {
	alpha<-0.1
}

plot(data$realInteractionDistance,data$observedSignal,ylim=c(allData.y.lim.bottom,allData.y.lim.top),main=paste(inputFileName,"1% top/bottom trim","C Scatter Plot - All Distances",sep="\n"),xlab="Genomic Distance (bp)",ylab="C counts",type="n")
points(data$realInteractionDistance,data$observedSignal,col=rgb(0,0,0,alpha))
lines(data$realInteractionDistance,data$loessExpectedValue,col="red",lwd=3)
lines(data$realInteractionDistance,data$loessExpectedValue+data$loessExpectedStdev,col="red",lwd=1,lty=2) # plot data$loessExpectedValue + data$loessExpectedStdev
lines(data$realInteractionDistance,data$loessExpectedValue-data$loessExpectedStdev,col="red",lwd=1,lty=2) # plot data$loessExpectedValue - data$loessExpectedStdev

legend("topright", legend = c("loess weighted average", "loess weighted stdev"),lty=1:2,lwd=3:1,xjust=1,col=c("red","red"),yjust=1)

# second plot - zoom in
allData.y.topIndex<-floor(allData.y.size*0.975)-1
allData.y.bottomIndex<-ceiling(allData.y.size*0.025)+1
allData.y.lim.top<-allData.y[allData.y.topIndex]
allData.y.lim.bottom<-allData.y[allData.y.bottomIndex]

plot(data$realInteractionDistance,data$observedSignal,ylim=c(allData.y.lim.bottom,allData.y.lim.top),main=paste(inputFileName,"5% top/bottom trim","C Scatter Plot - All Distances",sep="\n"),xlab="Genomic Distance (bp)",ylab="C counts",type="n")
points(data$realInteractionDistance,data$observedSignal,col=rgb(0,0,0,alpha))
lines(data$realInteractionDistance,data$loessExpectedValue,col="red",lwd=3)
lines(data$realInteractionDistance,data$loessExpectedValue+data$loessExpectedStdev,col="red",lwd=1,lty=2) # plot data$loessExpectedValue + data$loessExpectedStdev
lines(data$realInteractionDistance,data$loessExpectedValue-data$loessExpectedStdev,col="red",lwd=1,lty=2) # plot data$loessExpectedValue - data$loessExpectedStdev
legend("topright", legend = c("loess weighted average", "loess weighted stdev"),lty=1:2,lwd=3:1,xjust=1,col=c("red","red"),yjust=1)

dev.off()