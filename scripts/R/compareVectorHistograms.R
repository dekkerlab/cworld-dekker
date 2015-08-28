args <- commandArgs(TRUE)
wd<-args[1]
vector1File<-args[2]
vector2File<-args[3]
sampleName<-args[4]

setwd(wd)

vector1.name<-basename(vector1File)
vector1.data<-read.table(vector1File,header=T,sep="\t")
vector1.data.x<-as.numeric(vector1.data[,1])

vector2.name<-basename(vector2File)
vector2.data<-read.table(vector2File,header=T,sep="\t")
vector2.data.x<-as.numeric(vector2.data[,1])

tmp.data=c(vector1.data.x,vector2.data.x)
tmp.data<-sort(tmp.data)
tmp.data.size<-length(tmp.data)
tmp.data.topIndex=floor(tmp.data.size*.999)
tmp.data.bottomIndex=ceiling(tmp.data.size*.001)
tmp.data.max<-tmp.data[tmp.data.topIndex]
tmp.data.min<-tmp.data[tmp.data.bottomIndex]

if(tmp.data.min == tmp.data.max) {
	tmp.data.max<-tmp.data.min+1
}

tmp.data.range<-(tmp.data.max-tmp.data.min)
tmp.data.binWidth<-(tmp.data.range/100)

vector1.data.trimmed<-vector1.data.x[vector1.data.x>=tmp.data.min & vector1.data.x<=tmp.data.max]
vector2.data.trimmed<-vector2.data.x[vector2.data.x>=tmp.data.min & vector2.data.x<=tmp.data.max]

bins<-seq(tmp.data.min,tmp.data.max,by=tmp.data.binWidth)
bins.rounded<-seq(round(tmp.data.min,digits=1),round(tmp.data.max,digits=1),by=round(tmp.data.binWidth,digits=1))

pngfile<-paste(sampleName,".png",sep='')
png(pngfile,height=600,width=800)

hist(vector1.data.trimmed,freq=FALSE,breaks=bins,col=rgb(1,0,0,0.4),main=paste("Histogram",sampleName,sep="\n"),right=FALSE,xaxt="n")
hist(vector2.data.trimmed,freq=FALSE,add=T,breaks=bins,col=rgb(0,1,0,0.4),right=FALSE,xaxt="n")
abline(v=0,lwd=2,lty=1,col="black")
legend("topright", legend = c(vector1.name,vector2.name),pch=16,col=c("red","green"),xjust=1,yjust=1)

axis(1,bins.rounded)
dev.off()
