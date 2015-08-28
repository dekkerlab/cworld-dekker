args <- commandArgs(TRUE)
inputFile<-args[1]
signalColumn<-as.numeric(args[2])
nBreaks<-as.numeric(args[3])
skip<-as.numeric(args[4])
minX<-as.numeric(args[5])
maxX<-as.numeric(args[6])

if(is.na(skip)) {
	skip<-0
}

inputFileName<-basename(inputFile)
myData<-read.table(inputFile,skip=skip,header=F,sep="\t",comment.char="#")
nCols<-length(myData)

if(is.na(signalColumn)) {
	signalColumn<-1
} else {
	if(signalColumn > nCols) {
		signalColumn<-1
	}
}

myData<-subset(myData,myData[,signalColumn]!="NA" & myData[,signalColumn]!="NaN")
x<-as.numeric(myData[,signalColumn])

x.tmp<-sort(x)
x.tmp.size<-length(x.tmp)
x.tmp.topIndex=floor(x.tmp.size*0.95)
x.tmp.bottomIndex=ceiling(x.tmp.size*0.05)+1
x.tmp.max<-x.tmp[x.tmp.topIndex]
x.tmp.min<-x.tmp[x.tmp.bottomIndex]
rm(x.tmp)

inputFileName<-basename(inputFile)

pngfile<-paste(inputFileName,".png",sep='')
png(pngfile,height=600,width=800)

par(mar=c(5,5,10,5))

meanX<-signif(mean(x),8)
sdX<-signif(sd(x),8)
medianX<-signif(median(x),8)
iqrX<-signif(IQR(x),8)
lowerIQR<-signif((medianX-(1.5*iqrX)),8)
upperIQR<-signif((medianX+(1.5*iqrX)),8)

if(is.na(minX)) {
	minX<-x.tmp.min
}

if(is.na(maxX)) {
	maxX<-x.tmp.max
}

x<-x[x>=minX]
x<-x[x<=maxX]

if(is.nan(nBreaks)) {
	nBreaks<-10
	if(length(x) > 100) { nBreaks<-ceiling(length(x)/25) }
}
xBound<-max(abs(minX),abs(maxX))

breaks<-seq(from=minX,to=maxX,length.out=nBreaks)

label<-paste("nBreaks=",nBreaks," mean=",meanX," sdX=",sdX,"\n","median=",medianX," IQR=",iqrX,"\n","lowerIQR=",lowerIQR," upperIQR=",upperIQR,sep="")

hist(x,freq=FALSE,breaks=breaks,col=rgb(0.25,0.25,0.25),border=rgb(0.25,0.25,0.25),prob=TRUE,main=label,right=FALSE,xlim=c(minX,maxX))
res<-hist(x,breaks=breaks,col=rgb(0.25,0.25,0.25),border=rgb(0.25,0.25,0.25),prob=TRUE,main=label,right=FALSE,xlim=c(minX,maxX))

#cbind(res$mids,res$counts)

lines(density(x), col="blue", lwd=2)
lines(density(x, adjust=2), lty="dotted", col="green", lwd=2) 
abline(v=minX,lwd=2,lty=2,col="red")
abline(v=maxX,lwd=2,lty=2,col="red")

dev.off()

