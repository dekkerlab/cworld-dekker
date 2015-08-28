args <- commandArgs(TRUE)
wd<-args[1]
replicateFile<-args[2]
sampleFile<-args[3]
sampleName<-args[4]

setwd(wd)

rep.data<-read.table(replicateFile,header=T,sep="\t")
rep.data.x<-as.numeric(rep.data[,1])

sample.data<-read.table(sampleFile,header=T,sep="\t")
sample.data.x<-as.numeric(sample.data[,1])

tmp.data=c(rep.data.x,sample.data.x)
tmp.data<-sort(tmp.data)
tmp.data.size<-length(tmp.data)
tmp.data.topIndex=floor(tmp.data.size*.999)
tmp.data.bottomIndex=floor(tmp.data.size*.001)+1
tmp.data.max<-tmp.data[tmp.data.topIndex]
tmp.data.min<-tmp.data[tmp.data.bottomIndex]

if(tmp.data.min == tmp.data.max) {
	tmp.data.max<-tmp.data.min+1
}

tmp.data.range<-(tmp.data.max-tmp.data.min)
tmp.data.binWidth<-(tmp.data.range/100)

rep.data.trimmed<-rep.data.x[rep.data.x>=tmp.data.min & rep.data.x<=tmp.data.max]
sample.data.trimmed<-sample.data.x[sample.data.x>=tmp.data.min & sample.data.x<=tmp.data.max]

bins<-seq(tmp.data.min,tmp.data.max,by=tmp.data.binWidth)
bins.rounded<-seq(round(tmp.data.min,digits=1),round(tmp.data.max,digits=1),by=round(tmp.data.binWidth,digits=1))

pngfile<-paste(sampleName,".png",sep='')
png(pngfile,height=600,width=800)

hist(rep.data.trimmed,freq=FALSE,breaks=bins,col=rgb(1,0,0,0.4),main=paste("Score Histogram",sampleName,sep="\n"),right=FALSE,xaxt="n")
hist(sample.data.trimmed,freq=FALSE,add=T,breaks=bins,col=rgb(0,1,0,0.4),right=FALSE,xaxt="n")
abline(v=0,lwd=2,lty=1,col="black")
legend("topright", legend = c("replicate differences","sample differences"),pch=16,col=c("red","green"),xjust=1,yjust=1)

axis(1,bins.rounded)
dev.off()
