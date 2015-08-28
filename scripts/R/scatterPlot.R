
args <- commandArgs(TRUE)
inputFile_1<-args[1]
inputFile_2<-args[2]
name<-args[3]
dataColumn<-as.numeric(args[4])
skip<-as.numeric(args[5])
outlierFraction<-as.numeric(args[6])
correlationMode<-args[7]

if((outlierFraction > 1) | (outlierFraction < 0)) {
	outlierFraction<-0.05
}

name_1<-basename(inputFile_1)
name_2<-basename(inputFile_2)

if(length(name_1) > 60) {
	name_1<-substr(name_1, 1, 60)
}
if(length(name_2) > 60) {
	name_2<-substr(name_2, 1, 60)
}

plotFile<-paste(name,".",correlationMode,sep="")
pngfile<-paste(plotFile,".png",sep='')
png(pngfile,height=600,width=600)

tmpData_1<-read.table(inputFile_1,skip=skip,header=F,sep="\t",stringsAsFactors=FALSE,comment.char="#")
dataVector_1=tmpData_1[,eval(dataColumn)]

tmpData_2<-read.table(inputFile_2,skip=skip,header=F,sep="\t",stringsAsFactors=FALSE,comment.char="#")
dataVector_2=tmpData_2[,eval(dataColumn)]	

scatter<-cbind(dataVector_1,dataVector_2)
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
r<-rmatrix[2]
fit<-lm(y~x,data=scatter[good,])

plot(scatter, type = "n",main=paste(name,"\n",correlationMode,"\noutlier=",outlierFraction," & ","r=",r),cex=0.5,xlab=name_1,ylab=name_2)
points(scatter[-good,], col = "black", pch = 21, bg = "black", cex = 0.8)
points(scatter[good,], col = "red", pch = 21, bg = "red", cex = 0.8)
abline(fit, col = "blue", lwd = 2)

garbage <- dev.off()

cat(paste(name,"\t",r,"\n",sep=""))
