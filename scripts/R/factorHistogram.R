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
x.binWidth<-((x.bound*2)/100)
x.bins<-seq(-x.bound,x.bound,by=x.binWidth)

# calculate trimming
x.tmp<-sort(x)
x.tmp.size<-length(x.tmp)
x.tmp.topIndex=floor(x.tmp.size*.75)
x.tmp.bottomIndex=floor(x.tmp.size*.25)+1
x.tmp.min<-floor(x.tmp[x.tmp.bottomIndex])
x.tmp.max<-ceiling(x.tmp[x.tmp.topIndex])
# build the trimmed array
x.trimmed<-x.tmp[x.tmp>=x.tmp.min & x.tmp<=x.tmp.max]
x.trimmed.min<-floor(min(x.trimmed))
x.trimmed.max<-ceiling(max(x.trimmed))
x.trimmed.bound<-max(abs(x.trimmed.min),abs(x.trimmed.max))
if(max(abs(min(x.trimmed)),abs(max(x.trimmed))) < 0.5) { 
	x.trimmed.bound<-0.5
}
if(max(abs(min(x.trimmed)),abs(max(x.trimmed))) < 0.1) { 
	x.trimmed.bound<-0.1
}
x.trimmed.binWidth<-((x.trimmed.bound*2)/100)
x.trimmed.bins<-seq(-x.trimmed.bound,x.trimmed.bound,by=x.trimmed.binWidth)

pngfile<-paste(inputFile,".png",sep='')
png(pngfile,height=600,width=1100)

par(mfrow=c(2,1))

hist(x,breaks=x.bins,col="grey",main=paste("Factor Histogram",inputFileName,sep="\n"),right=FALSE,xaxt="n",xlim=c(-x.bound,x.bound),xlab="factors",ylab="frequency")
abline(v=0,col="black",lty=2,lwd=2)
axis(1,seq(from=-x.bound,to=x.bound,by=((x.bound*2)/10)))

hist(x.trimmed,breaks=x.trimmed.bins,col="grey",main=paste("Factor Histogram",inputFileName,sep="\n"),right=FALSE,xaxt="n",xlim=c(-x.trimmed.bound,x.trimmed.bound),xlab="factors",ylab="frequency")
abline(v=0,col="black",lty=2,lwd=2)
axis(1,seq(from=-x.trimmed.bound,to=x.trimmed.bound,by=((x.trimmed.bound*2)/10)))

dev.off()