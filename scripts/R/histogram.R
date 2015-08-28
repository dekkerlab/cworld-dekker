args <- commandArgs(TRUE)
wd<-args[1]
inputFile<-args[2]

setwd(wd)

data<-read.table(inputFile,header=F,sep="\t")
data.x<-as.numeric(data[,1])

x.tmp<-sort(data.x)
x.tmp.size<-length(x.tmp)
x.tmp.topIndex=floor(x.tmp.size*.99)
x.tmp.bottomIndex=floor(x.tmp.size*.01)+1
x.tmp.max<-x.tmp[x.tmp.topIndex]
x.tmp.min<-x.tmp[x.tmp.bottomIndex]

if(x.tmp.min == x.tmp.max) {
	x.tmp.max<-x.tmp.min+1
}

x.tmp.range<-(x.tmp.max-x.tmp.min)
x.tmp.binWidth<-(x.tmp.range/100)

x.trimmed<-x.tmp[x.tmp>=x.tmp.min & x.tmp<=x.tmp.max]
bins<-seq(x.tmp.min,x.tmp.max,by=x.tmp.binWidth)

x.trimmed.neg<-x.trimmed[x.trimmed<0]
x.trimmed.pos<-x.trimmed[x.trimmed>0]

x.trimmed.neg.mean<-mean(x.trimmed.neg)
x.trimmed.neg.sd<-sd(x.trimmed.neg)
x.trimmed.neg.mean.minus.sd=x.trimmed.neg.mean-x.trimmed.neg.sd
x.trimmed.neg.mean.plus.sd=x.trimmed.neg.mean+x.trimmed.neg.sd


x.trimmed.pos.mean<-mean(x.trimmed.pos)
x.trimmed.pos.sd<-sd(x.trimmed.pos)
x.trimmed.pos.mean.minus.sd=x.trimmed.pos.mean-x.trimmed.pos.sd
x.trimmed.pos.mean.plus.sd=x.trimmed.pos.mean+x.trimmed.pos.sd

pngfile<-paste(inputFile,".png",sep='')
png(pngfile,height=400,width=600)
hist(x.trimmed,breaks=bins,col="grey",main=paste("Score Histogram",inputFile,sep="\n"),right=FALSE,xaxt="n")
abline(v=0,lwd=2,lty=2) # neg/pos seperator
abline(v=x.trimmed.neg.mean,lwd=2,col="red") 
abline(v=x.trimmed.neg.mean.minus.sd,lwd=1,lty=2,col="red")
abline(v=x.trimmed.neg.mean.plus.sd,lwd=1,lty=2,col="red")
abline(v=x.trimmed.pos.mean,lwd=2,col="blue")
abline(v=x.trimmed.pos.mean.minus.sd,lwd=1,lty=2,col="blue")
abline(v=x.trimmed.pos.mean.plus.sd,lwd=1,lty=2,col="blue")
axis(1,bins)
dev.off()

