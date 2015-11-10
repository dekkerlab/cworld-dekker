options(bitmapType='cairo')

args <- commandArgs(TRUE)
dir<-args[1]
inputFile<-args[2]
inputFileName<-args[3]

setwd(dir)

data<-read.table(inputFile,header=F,sep="\t")

x<-data[,2]

# calcuate for full data range
x.min<-floor(min(x))
x.max<-ceiling(max(x))
x.bound<-max(abs(x.min),abs(x.max))
if(x.bound < 1) {
	x.bound<-1
}
x.binWidth<-((x.bound*2)/300)
x.bins<-seq(-x.bound,x.bound,by=x.binWidth)

# calculate trimming
x.tmp<-sort(x)
x.tmp.size<-length(x.tmp)
x.tmp.topIndex=floor(x.tmp.size*.75)
x.tmp.bottomIndex=floor(x.tmp.size*.25)+1
x.tmp.min<-floor(x.tmp[x.tmp.bottomIndex])
x.tmp.max<-ceiling(x.tmp[x.tmp.topIndex])

pngfile<-paste(inputFile,".png",sep='')
png(pngfile,height=500,width=900)

hist(x,breaks=x.bins,col="grey",main=paste("Factor Histogram",inputFileName,sep="\n"),right=FALSE,xaxt="n",xlim=c(-x.bound,x.bound),xlab="factors",ylab="frequency")
abline(v=0,col="black",lty=2,lwd=2)
axis(1,seq(from=-x.bound,to=x.bound,by=((x.bound*2)/10)))

dev.off()