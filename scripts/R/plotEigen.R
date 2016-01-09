options(bitmapType='cairo')

plotEigen <- function(x,y,plotTitle) {
	
	ymin<-min(y,na.rm=TRUE)
	ymax<-max(y,na.rm=TRUE)
	
	xStart<-min(x)
	xEnd<-max(x)
	
	ybound<-max(abs(ymin),abs(ymax))
	if(ymin < 0) {
		ybounds<-c(-ybound,ybound)
	} else {
		ybounds<-c(0,ybound)
	}
	
	plot(x,y,type="n",main=plotTitle,ylab="eigen value",xlab="genomic coordinates (bp)",ylim=ybounds)

	for (i in 1:(length(x)-1) ) {
		tmpX<-x[i]
		tmpY<-y[i]
		nextX<-x[i+1]
		nextY<-y[i+1]
		if(!is.nan(tmpY) & !is.nan(nextY)) {
			if((tmpY*nextY) >= 0) { # same sign? +/+ or -/0
				polygon(c(tmpX,nextX,nextX,tmpX), c(tmpY,nextY,0,0), col = ifelse(nextY >= 0,'red','blue'), border = NA)
			} else { # if we are about to cross 0 line - then special case polygons
				polygon(c(tmpX,nextX,tmpX), c(tmpY,0,0), col = ifelse(tmpY >= 0,'red','blue'), border = NA)
				polygon(c(tmpX,nextX,tmpX), c(tmpY,0,0), col = ifelse(tmpY >= 0,'red','blue'), border = NA)
			}
		} else {
			polygon(c(tmpX,nextX,nextX,tmpX), c(-ybound,-ybound,ybound,ybound), col=rgb(0.75,0.75,0.75,0.5), border = NA)
		}
	}

	abline(h=0,lwd=1,lty=2,col="black")
}

plotGeneDensity <- function(x,y,eigen1,plotTitle) {
	ymin<-min(y,na.rm=TRUE)
	ymax<-max(y,na.rm=TRUE)
	ybound<-max(abs(ymin),abs(ymax))
	if(ymin < 0) {
		ybounds<-c(-ybound,ybound)
	} else {
		ybounds<-c(0,ybound)
	}
	    
    plot(x,y,type="n",main=plotTitle,ylab="gene density",xlab="genomic coordinates (bp)",ylim=ybounds)
	
	for (i in 1:(length(x)-1) ) {
		tmpX<-x[i]
		tmpY<-y[i]
        tmpEigen1<-eigen1[i]
		nextX<-x[i+1]
		nextY<-y[i+1]
		
		if(!is.nan(tmpY) & !is.nan(nextY)) {
			polygon(c(tmpX,nextX,nextX,tmpX), c(tmpY,nextY,0,0), col = ifelse(tmpEigen1 >= 0,'black','darkgray'), border = NA)
		} else {
			polygon(c(tmpX,nextX,nextX,tmpX), c(-ybound,-ybound,ybound,ybound), col=rgb(255,165,0,100,maxColorValue=255), border = NA)
		}
	}

	abline(h=0,lwd=1,lty=2,col="black")
}


args <- commandArgs(TRUE)
cwd<-args[1]
inputFile<-args[2]
jobName<-args[3]

setwd(cwd)

data<-read.table(inputFile,header=T,sep="\t",comment.char="")
x<-data$index

pngfile<-paste(inputFile,".png",sep='')
png(pngfile,height=1200,width=1200)

par(mfrow= c(4, 1))
par(cex=1.1)

plotGeneDensity(x,data$geneDensity,data$eigen1,paste(jobName,"geneDensity",sep="\n"))

plotEigen(x,data$eigen1,paste("prinicple component #1 - evr [",(mean(data$eigen1evr)*100),"%]",sep=""))
plotEigen(x,data$eigen2,paste("prinicple component #2 - evr [",(mean(data$eigen2evr)*100),"%]",sep=""))
plotEigen(x,data$eigen3,paste("prinicple component #3 - evr [",(mean(data$eigen3evr)*100),"%]",sep=""))

dev.off()

