options(bitmapType='cairo')

args <- commandArgs(TRUE)
dir<-args[1]
inputFile<-args[2]
headerSize<-as.numeric(args[3])
headerSpacing<-as.numeric(args[4])
binSize<-as.numeric(args[5])
binSizeDistance<-as.numeric(args[6])
imageWidth<-as.numeric(args[7])
imageHeight<-ceiling(imageWidth/8)
yBound<-as.numeric(args[8])
transparentBGFlag<-as.numeric(args[9])

setwd(dir)

insulationData<-read.table(inputFile,header=T,sep="\t")

xStart<-min(insulationData$start)
xEnd<-max(insulationData$end)
xRange<-(xEnd-xStart)
xBinStart<-min(insulationData$binStart)
xBinEnd<-max(insulationData$binEnd)
xBinRange<-(xBinEnd-xBinStart)

insulationData<-subset(insulationData,insulationData$insulationScore!="NA") # remove NAs from the scores
insulationData<-subset(insulationData,insulationData$delta!="NA") # remove NAs from the scores
insulationData<-subset(insulationData,insulationData$deltaSquare!="NA") # remove NAs from the scores

insulationData<-insulationData[with(insulationData, order(start)), ]

x<-insulationData$binMidpoint
y<-insulationData$insulationScore
delta<-insulationData$delta
deltaSquare<-insulationData$deltaSquare

yStart<-min(y)
yEnd<-max(y)

if(yBound == 0) {
	yBound<-ceiling(max(abs(yStart),abs(yEnd)))
	if(yBound < 1) {	
		yBound<-1
	}
}

yStart<--yBound
yEnd<-yBound
yRange<-(yEnd-yStart)

pos_yBound <- yBound
neg_yBound <- -yBound

# saturate plot at specified yBound(s)
y[y > pos_yBound] <- pos_yBound
y[y < neg_yBound] <- neg_yBound

inputFile<-gsub(".gz", "", inputFile)
insulationData.boundaries<-read.table(paste(inputFile,".boundaries",sep=""),header=T,sep="\t")

x.start.boundaries<-insulationData.boundaries$binStart
x.end.boundaries<-insulationData.boundaries$binEnd
y.boundaries<-insulationData.boundaries$boundaryStrength

# set minimum image height
if(imageHeight < 325) {
	imageHeight<-325
}

# debug image

# png file
pngfile<-paste(inputFile,".debug.png",sep='')
if(transparentBGFlag == 1) {
	png(pngfile,height=imageHeight,width=imageWidth,bg="transparent")
} else {
	png(pngfile,height=imageHeight,width=imageWidth)
}

par(mar=c(4, 4, 4, 4) + 0.1)

plot(x,y,main=paste(inputFile,"\n","binSize=",binSize,"|",binSizeDistance,sep=""),cex=0.5,col=rgb(0.25,0.25,0.25,0.5),xlab="genomic coordinates",ylab="insulation index",xlim=c(xBinStart,xBinEnd),ylim=c(-yBound,yBound),type="n",xaxt="n",yaxt="n")
axis(2,seq(from=-yBound,to=yBound,length.out=11))

if(length(x) > 1) {

	if(nrow(insulationData.boundaries) > 0) {
		for (i in 1:nrow(insulationData.boundaries)) {
			tmp.x.start.boundary<-x.start.boundaries[i]
			tmp.x.end.boundary<-x.end.boundaries[i]
			tmp.y.boundary<-y.boundaries[i]
			
			strength<-1
            
			#segments(tmp.x.boundary, yStart, tmp.x.boundary, yEnd, col=rgb(0,1,0,strength), lwd=1)
			rect(tmp.x.start.boundary, yStart, tmp.x.end.boundary, yEnd, col=rgb(0,1,0,strength), border=NA)
		}
	}
	
	for (i in 1:(length(x)-1) ) {
		tmpX<-x[i]
		tmpY<-y[i]
		nextX<-x[i+1]
		nextY<-y[i+1]
		
		if((i == 1) && ((tmpX-xBinStart) != 0)) {
			rect(xBinStart, yStart, tmpX, yEnd, col=rgb(0.75,0.75,0.75,0.5),border=FALSE,lwd=0)
		}
		
		if((nextX-tmpX) == 1) {
			segments(tmpX,tmpY,nextX,nextY,col="black",lwd=1)
		} else {
			rect(tmpX, yStart, nextX, yEnd, col=rgb(0.75,0.75,0.75,0.5),border=FALSE,lwd=1)
		}
		
		if((i == (length(x)-1)) && ((xBinEnd-tmpX) != 0)) {
			rect(tmpX, yStart, xBinEnd, yEnd, col=rgb(0.75,0.75,0.75,0.5),border=FALSE,lwd=1)
		}
	}
	
	# draw saturation blobs
	for (i in 1:(length(x)-1) ) {
		tmpX<-x[i]
		tmpY<-y[i]
		nextX<-x[i+1]
		nextY<-y[i+1]
		
		if((nextX-tmpX) == 1) {
			if(nextY >= yBound) {
				segments(tmpX,yBound,nextX,yBound,col="purple",lwd=4)
			} else if(nextY <= -yBound) {
				segments(tmpX,-yBound,nextX,-yBound,col="purple",lwd=4)
			}
		}
		
	}

	abline(h=0,lwd=1,lty=2,col="black")
	abline(v=xBinStart,col="red",lwd=2,lty=2)
	abline(v=xBinEnd,col="red",lwd=2,lty=2)
	
}

par(new=T)

deltaStart<-min(delta)
deltaEnd<-max(delta)
deltaBound<-max(abs(deltaStart),abs(deltaEnd))

deltaStart<--deltaBound
deltaEnd<-deltaBound
deltaRange<-(deltaEnd-deltaStart)

plot(x,delta,xlim=c(xBinStart,xBinEnd),ylim=c(-deltaBound,deltaBound),xlab="",ylab="",axes=FALSE,xaxt="n",yaxt="n",type="n")
axis(4,seq(from=-deltaBound,to=deltaBound,length.out=11))

if(length(x) > 1) {
	for (i in 1:(length(x)-1) ) {
		tmpX<-x[i]
		tmpY<-delta[i]
		nextX<-x[i+1]
		nextY<-delta[i+1]
		
		if((nextX-tmpX) == 1) {
			segments(tmpX,tmpY,nextX,nextY,col="blue",lwd=1)
		}
	}
}

mtext("delta",side=4)

par(new=T)
plot(x,deltaSquare,xlim=c(xBinStart,xBinEnd),ylim=c(-2,2),xlab="",ylab="",axes=FALSE,xaxt="n",yaxt="n",type="n")
lines(x,deltaSquare,col="red",lwd=1)

axis(1,at=seq(from=xBinStart,to=xBinEnd,length.out=11),labels=seq(from=xStart,to=xEnd,length.out=11))

dev.off()

# production image
# png file
pngfile<-paste(inputFile,".png",sep='')
if(transparentBGFlag == 1) {
	png(pngfile,height=imageHeight,width=imageWidth,bg="transparent")
} else {
	png(pngfile,height=imageHeight,width=imageWidth)
}

par(mar=c(4, 4, 4, 4) + 0.1)

plot(x,y,main=paste(inputFile,"\n","binSize=",binSize,"|",binSizeDistance,sep=""),cex=0.5,col=rgb(0.25,0.25,0.25,0.5),xlab="genomic coordinates",ylab="insulation index",xlim=c(xBinStart,xBinEnd),ylim=c(-yBound,yBound),type="n",xaxt="n",yaxt="n")
axis(2,seq(from=-yBound,to=yBound,length.out=11))

if(length(x) > 1) {

	if(nrow(insulationData.boundaries) > 0) {
		for (i in 1:nrow(insulationData.boundaries)) {
			tmp.x.start.boundary<-x.start.boundaries[i]
			tmp.x.end.boundary<-x.end.boundaries[i]
			tmp.y.boundary<-y.boundaries[i]
			
			strength<-1
			#segments(tmp.x.boundary, yStart, tmp.x.boundary, yEnd, col=rgb(0,1,0,strength), lwd=1)
			rect(tmp.x.start.boundary, yStart, tmp.x.end.boundary, yEnd, col=rgb(0,1,0,strength), border=NA)
		}
	}
	
	for (i in 1:(length(x)-1) ) {
		tmpX<-x[i]
		tmpY<-y[i]
		nextX<-x[i+1]
		nextY<-y[i+1]
		
		if((i == 1) && ((tmpX-xBinStart) != 0)) {
			rect(xBinStart, yStart, tmpX, yEnd, col=rgb(0.75,0.75,0.75,0.5),border=FALSE,lwd=0)
		}
		
		if((nextX-tmpX) == 1) {
			segments(tmpX,tmpY,nextX,nextY,col="black",lwd=1)
		} else {
			rect(tmpX, yStart, nextX, yEnd, col=rgb(0.75,0.75,0.75,0.5),border=FALSE,lwd=1)
		}
		
		if((i == (length(x)-1)) && ((xBinEnd-tmpX) != 0)) {
			rect(tmpX, yStart, xBinEnd, yEnd, col=rgb(0.75,0.75,0.75,0.5),border=FALSE,lwd=1)
		}
	}
	
	# draw saturation blobs
	for (i in 1:(length(x)-1) ) {
		tmpX<-x[i]
		tmpY<-y[i]
		nextX<-x[i+1]
		nextY<-y[i+1]
		
		if((nextX-tmpX) == 1) {
			if(nextY >= yBound) {
				segments(tmpX,yBound,nextX,yBound,col="purple",lwd=4)
			} else if(nextY <= -yBound) {
				segments(tmpX,-yBound,nextX,-yBound,col="purple",lwd=4)
			}
		}
		
	}

	abline(h=0,lwd=1,lty=2,col="black")
	abline(v=xBinStart,col="red",lwd=2,lty=2)
	abline(v=xBinEnd,col="red",lwd=2,lty=2)
	
}

axis(1,at=seq(from=xBinStart,to=xBinEnd,length.out=11),labels=seq(from=xStart,to=xEnd,length.out=11))

dev.off()





# pdf file
pdffile<-paste(inputFile,".pdf",sep='')
pdf(pdffile)

par(mar=c(4, 4, 4, 4) + 0.1)

plot(x,y,main=paste(inputFile,"\n","binSize=",binSize,"|",binSizeDistance,sep=""),cex=0.5,col=rgb(0.25,0.25,0.25,0.5),xlab="genomic coordinates",ylab="insulation index",xlim=c(xBinStart,xBinEnd),ylim=c(-yBound,yBound),type="n",xaxt="n",yaxt="n")
axis(2,seq(from=-yBound,to=yBound,length.out=11))

if(length(x) > 1) {

	if(nrow(insulationData.boundaries) > 0) {
		for (i in 1:nrow(insulationData.boundaries)) {
			tmp.x.start.boundary<-x.start.boundaries[i]
			tmp.x.end.boundary<-x.end.boundaries[i]
			tmp.y.boundary<-y.boundaries[i]
			
			strength<-1
			#segments(tmp.x.boundary, yStart, tmp.x.boundary, yEnd, col=rgb(0,1,0,strength), lwd=1)
			rect(tmp.x.start.boundary, yStart, tmp.x.end.boundary, yEnd, col=rgb(0,1,0,strength), border=NA)
		}
	}
	
	for (i in 1:(length(x)-1) ) {
		tmpX<-x[i]
		tmpY<-y[i]
		nextX<-x[i+1]
		nextY<-y[i+1]
		
		if((i == 1) && ((tmpX-xBinStart) != 0)) {
			rect(xBinStart, yStart, tmpX, yEnd, col=rgb(0.75,0.75,0.75,0.5),border=FALSE,lwd=0)
		}
		
		if((nextX-tmpX) == 1) {
			segments(tmpX,tmpY,nextX,nextY,col="black",lwd=1)
		} else {
			rect(tmpX, yStart, nextX, yEnd, col=rgb(0.75,0.75,0.75,0.5),border=FALSE,lwd=1)
		}
		
		if((i == (length(x)-1)) && ((xBinEnd-tmpX) != 0)) {
			rect(tmpX, yStart, xBinEnd, yEnd, col=rgb(0.75,0.75,0.75,0.5),border=FALSE,lwd=1)
		}
	}
	
	# draw saturation blobs
	for (i in 1:(length(x)-1) ) {
		tmpX<-x[i]
		tmpY<-y[i]
		nextX<-x[i+1]
		nextY<-y[i+1]
		
		if((nextX-tmpX) == 1) {
			if(nextY >= yBound) {
				segments(tmpX,yBound,nextX,yBound,col="purple",lwd=4)
			} else if(nextY <= -yBound) {
				segments(tmpX,-yBound,nextX,-yBound,col="purple",lwd=4)
			}
		}
		
	}

	abline(h=0,lwd=1,lty=2,col="black")
	abline(v=xBinStart,col="red",lwd=2,lty=2)
	abline(v=xBinEnd,col="red",lwd=2,lty=2)
	
}

axis(1,at=seq(from=xBinStart,to=xBinEnd,length.out=11),labels=seq(from=xStart,to=xEnd,length.out=11))

dev.off()