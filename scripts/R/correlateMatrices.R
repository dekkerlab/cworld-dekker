options(bitmapType='cairo')

args <- commandArgs(TRUE)
dir<-args[1]
inputFile<-args[2]
name_1<-args[3]
name_2<-args[4]
correlationMode<-args[5]
outlierFraction<-as.numeric(args[6])
ymin<-as.numeric(args[7])
ymax<-as.numeric(args[8])
xmin<-as.numeric(args[9])
xmax<-as.numeric(args[10])

if(is.na(outlierFraction) | (outlierFraction > 1) | (outlierFraction < 0)) {
	outlierFraction<-0.05
}

setwd(dir)

plotFile<-inputFile
plotFile<-basename(plotFile)
plotFile<-gsub(".gz", "", plotFile)
plotFile<-paste(plotFile,".",correlationMode,sep="")

name=paste(name_1,"\n",name_2,sep="")
pngfile<-paste(plotFile,".png",sep='')
png(pngfile,height=600,width=600)

if((file.info(inputFile)[[1]][1]) > 0) {
	
	myData<-read.table(inputFile,header=TRUE,sep="\t")
    origRows=nrow(myData)
    
    scatter<-cbind(myData$cScore_1,myData$cScore_2)
    scatter<-subset(scatter,scatter[,1]!="NA")
    scatter<-subset(scatter,scatter[,1]!="nan")
    scatter<-subset(scatter,scatter[,1]!="NaN")
    scatter<-subset(scatter,scatter[,2]!="NA")
    scatter<-subset(scatter,scatter[,2]!="nan")
    scatter<-subset(scatter,scatter[,2]!="NaN")

    scatter <- data.frame(scatter)
    names(scatter) <- c("x","y")
    
    res <- resid(mod <- lm(y ~ x,data=scatter))
    res.qt <- quantile(res, probs = c(outlierFraction,(1-outlierFraction)))
    good <- which(res >= res.qt[1] & res <= res.qt[2])
    rmatrix<-cor(scatter[good,],method=correlationMode)
    r<-round(rmatrix[2],6)
    fit<-lm(y~x,data=scatter[good,])

    if(is.na(ymin)) { ymin<-min(scatter[good,]$y) }
    if(is.na(ymax)) { ymax<-max(scatter[good,]$y) }
    if(is.na(xmin)) { xmin<-min(scatter[good,]$x) }
    if(is.na(xmax)) { xmax<-max(scatter[good,]$x) }
    
    plot(scatter, type = "n",main=paste(name,"\n",correlationMode," & outlier=",outlierFraction," & ","r=",r),cex=0.5,xlab=name_1,ylab=name_2,ylim=c(ymin,ymax),xlim=c(xmin,xmax))
    points(scatter[-good,], col = "black", pch = 21, bg = "black", cex = 0.8)
    points(scatter[good,], col = "red", pch = 21, bg = "red", cex = 0.8)
    abline(fit, col = "blue", lwd = 2)

} else {
	
	plot(0,0,cex=0.5,xlab=name_1,ylab=name_2,main="ERROR - Cannot perform correlation")

}

dev.off()