options(bitmapType='cairo')

args <- commandArgs(TRUE)
name<-args[1]
dataColumn<-as.numeric(args[2])
skip<-as.numeric(args[3])
minX<-as.numeric(args[4])
maxX<-as.numeric(args[5])

fileArray=args[6:length(args)]
nFiles=length(fileArray)

cat("found ",nFiles," files.\n")
for (i in 1:nFiles) {
	arg<-fileArray[i]
	cat("\t",i,"\t",arg,"\n")
}

allX <- vector()
myData <- list()
nameData <- list()
minSize<-0
for (i in 1:nFiles) {
	tmpFile<-fileArray[i]
	
	tmpData<-read.table(tmpFile,header=F,skip=skip,sep="\t",stringsAsFactors=FALSE,comment.char="#")
	dataVector=tmpData[,eval(dataColumn)]
	dataVector<-subset(dataVector,dataVector!="NaN")
	dataVector<-subset(dataVector,dataVector!="NA")
	dataVector<-subset(dataVector,dataVector!=".")
	
	dataVector<-as.numeric(dataVector)
	
	dataVector<-sort(dataVector)
	
	allX<-c(allX,dataVector)
	myData[[i]]<-dataVector
	
	if(minSize == 0) { minSize<-length(dataVector) }
	minSize<-min(minSize,length(dataVector))
	
	tmpName<-basename(tmpFile)
	if(nchar(tmpName) > 40) {
		tmpName<-substr(tmpName,0,40)
	}
	nameData[i]<-paste("(",length(dataVector),") ",tmpName,sep="")
	
	
}

allX.length<-length(allX)
allX.sorted<-sort(allX)

x.tmp.min<-floor(allX.sorted[ceiling(allX.length*0.001)])
x.tmp.max<-ceiling(allX.sorted[floor(allX.length*0.999)])

if(is.na(minX)) {
	minX<-x.tmp.min
}
if(is.na(maxX)) {
	maxX<-x.tmp.max
}

# png file

pngfile<-paste(name,".png",sep='')
png(pngfile,height=600,width=900)

par(mar=c(5,5,5,5))

lineColors<-rainbow(nFiles,alpha=0.05)

nBreaks<-10
if(minSize > 100) { nBreaks<-ceiling(minSize/25) }
hist(myData[[1]],breaks=nBreaks,freq=FALSE,right=FALSE,col=lineColors[1])

for (i in 1:nFiles) {
	x<-myData[[i]]
	tmpColor=lineColors[i]
	
	hist(x,breaks=nBreaks,freq=FALSE,add=T,right=FALSE,col=tmpColor)
}

dev.off()

