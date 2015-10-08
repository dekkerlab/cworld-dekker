options(bitmapType='cairo')

args <- commandArgs(TRUE)
wd<-args[1]
name<-args[2]

setwd(wd)

inputfile<-paste(name,sep='')
data<-read.table(inputfile,header=F,sep="\t")

dist<-data[,3]
raw<-data[,4]
loess<-data[,5]
stdev<-data[,5]+data[,6]

y<-raw
y2<-y
y2<-sort(y2)
ysize<-length(y2)
ysize<-floor(ysize*.995)
ysize<-y2[ysize]
ymin<-y2[1]

pngfile<-paste(name,".png",sep='')
png(pngfile,height=720,width=1280)
plot(dist,raw,ylim=c(0,ysize),main=paste(name,"5C Scatter Plot - All Distances",sep="\n"),xlab="Genomic Distance (bp)",ylab="5C counts")
lines(dist,loess,col="red",lwd=3)
lines(dist,stdev,col="red",lwd=1,lty=2)
legend("topright", legend = c("loess weighted average", "loess weighted stdev"),lty=1:2,lwd=3:1,xjust=1,col=c("red","red"),yjust=1)
dev.off()

