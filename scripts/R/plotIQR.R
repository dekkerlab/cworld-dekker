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

barData <- vector()
boxData <- list()
nameData <- list()
for (i in 1:nFiles) {
	tmpFile<-fileArray[i]
	
	tmpData<-read.table(tmpFile,skip=skip,header=F,sep="\t",stringsAsFactors=FALSE,comment.char="#")
	dataVector=vector()
	dataVector=tmpData[,eval(dataColumn)]
	dataVector<-subset(dataVector,dataVector!="NaN")
	dataVector<-subset(dataVector,dataVector!="NA")
	dataVector<-subset(dataVector,dataVector!=".")
	
	dataVector<-as.numeric(dataVector)
	
	tmpName<-basename(tmpFile)
	tmpName.arr<-unlist(strsplit(tmpName,"\\."))
	tmpName<-tmpName.arr[1]
	if(nchar(tmpName) > 50) {
		tmpName<-substr(tmpName,0,50)
	}
	
    # do some outlier trimming (1.5 IQR)
	tmpIQR<-IQR(dataVector)
	tmpMedian<-median(dataVector)
	dataVector <- sort(dataVector)
	nData<-length(dataVector)
	q1<-dataVector[floor(nData*0.25)]
	q3<-dataVector[floor(nData*0.75)]
	tmpTop<-tmpMedian+(1.5*tmpIQR)
	tmpBottom<-tmpMedian-(1.5*tmpIQR)
	
	cat(paste(tmpName,tmpIQR,"\n"))
	barData<-c(barData,tmpIQR)
	nameData<-c(nameData,paste("(",length(dataVector),") ",tmpName,sep=""))
	
	# do some outlier trimming (1.5 IQR)
	#dataVector<-subset(dataVector,dataVector>tmpBottom)
	#dataVector<-subset(dataVector,dataVector<tmpTop)
	
	#tmpIQR<-IQR(dataVector)
	
	#cat(paste(tmpName,tmpIQR,"\n"))
	#barData<-c(barData,tmpIQR)
	#nameData<-c(nameData,paste("(",length(dataVector),") ",tmpName,sep=""))
	#nameData[i]<-paste(tmpName,sep="")
	
}

height=325
width=800

# png files

pngfile<-paste(name,".png",sep='')
png(pngfile,height=325,width=800)

par(mar=c(15,5,2,2))
barplot(barData,names.arg=nameData,col=rainbow(length(barData)),ylim=c(minY,maxY),las=2)

dev.off()


# pdf files

pdffile<-paste(name,".pdf",sep='')
pdf(pdffile)

par(mar=c(15,5,2,2))
barplot(barData,names.arg=nameData,col="orange",ylim=c(minY,maxY),las=2)

dev.off()