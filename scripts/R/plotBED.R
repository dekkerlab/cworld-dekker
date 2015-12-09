options(bitmapType='cairo')

args <- commandArgs(TRUE)
cwd<-args[1]
inputFile<-args[2]
jobName<-args[3]
yLimit<-as.numeric(args[4])

setwd(cwd)

data<-read.table(inputFile,header=F,skip=1,sep="\t")

x<-data[,2]
y<-data[,5]


xStart<-min(x)
xEnd<-max(x)

ymin<-min(y,na.rm=TRUE)
ymax<-max(y,na.rm=TRUE)
ybound<-max(abs(ymin),abs(ymax))
if(!is.na(yLimit) && (yLimit != "NA")) {
	ybound<-yLimit
}

pngfile<-paste(inputFile,".png",sep='')

png(pngfile,height=225,width=1200)

plot(x,y,type="n",main=paste(jobName,sep=""),ylab="BED signal",xlab="genomic coordinates (bp)",xlim=c(xStart,xEnd),ylim=c(-ybound,ybound))

lines(x,y,col="blue",lwd=2)
abline(h=0,lwd=1,lty=2,col="black")

dev.off()