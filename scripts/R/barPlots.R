options(bitmapType='cairo')

args <- commandArgs(TRUE)
name<-args[1]
dataColumn<-as.numeric(args[2])
skip<-as.numeric(args[3])
minY<-as.numeric(args[4])
maxY<-as.numeric(args[5])

fileArray=args[6:length(args)]
nFiles=length(fileArray)

cat("found ",nFiles," files.\n")
for (i in 1:nFiles) {
	arg<-fileArray[i]
	cat("\t",i,"\t",arg,"\n")
}

allY <- vector()
boxData <- list()
nameData <- list()
for (i in 1:nFiles) {
	tmpFile<-fileArray[i]
	
	tmpData<-read.table(tmpFile,skip=skip,header=F,sep="\t",stringsAsFactors=FALSE,comment.char="#")
	dataVector=tmpData[,eval(dataColumn)]
	dataVector<-as.numeric(dataVector)
	
	dataVector<-subset(dataVector,dataVector!="NaN")
	dataVector<-subset(dataVector,dataVector!="NA")
	dataVector<-subset(dataVector,dataVector!=".")
	#dataVector<-abs(dataVector)
	
	allY<-c(allY,dataVector)
	boxData<-c(boxData,mean(dataVector))
	
	tmpName<-basename(tmpFile)
	tmpName.arr<-unlist(strsplit(tmpName,"\\."))
	tmpName<-tmpName.arr[1]
	print(tmpName)
	
	if(nchar(tmpName) > 75) {
		tmpName<-substr(tmpName,0,75)
	}
	nameData[i]<-paste("(",length(dataVector),") ",tmpName,sep="")
	nameData[i]<-paste(tmpName,sep="")
	
}

allY.length<-length(allY)
allY.sorted<-sort(allY)

y.tmp.min<-floor(allY.sorted[ceiling(allY.length*0.001)])
y.tmp.max<-ceiling(allY.sorted[floor(allY.length*0.999)])

if(is.na(minY)) {
	minY<-y.tmp.min
}
if(is.na(maxY)) {
	maxY<-y.tmp.max
}

#boxWidths<-rep((700/length(boxData)), length(boxData))
minBoxWidth<-20
boxWidths<-rep(minBoxWidth, length(boxData))
width<-(minBoxWidth*length(boxData))+400
if(width < 800) {
	width=800
	boxWidths<-rep(floor(400/length(boxData)),length(boxData))
}
height<-ceiling(width/3)
if(height < 600) {
	height<-600
}

# png files

pngfile<-paste(name,".png",sep='')
png(pngfile,height=height,width=width)

par(mar=c(25,5,2,2))
barplot(unlist(boxData),names=nameData,width=boxWidths,col=rainbow(length(boxData)),ylim=c(minY,maxY),las=2,yaxt="n")
axis(2,seq(from=minY,to=maxY,length.out=5))
abline(h=0,lty=2,lwd=1,col="black")

dev.off()

# pdf files

pdffile<-paste(name,".pdf",sep='')
pdf(pdffile)

par(mar=c(25,5,2,2))
barplot(unlist(boxData),names=nameData,width=boxWidths,col="orange",ylim=c(minY,maxY),las=2,yaxt="n")
axis(2,seq(from=minY,to=maxY,length.out=5))
abline(h=0,lty=2,lwd=1,col="black")

dev.off()
