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
	dataVector=tmpData[,eval(dataColumn)]
	dataVector<-subset(dataVector,dataVector!="NaN")
	dataVector<-subset(dataVector,dataVector!="NA")
	dataVector<-subset(dataVector,dataVector!=".")
	
	dataVector<-as.numeric(dataVector)

	tmpName<-basename(tmpFile)
	tmpName.arr<-unlist(strsplit(tmpName,"\\."))
	print(tmpName)
	tmpName<-tmpName.arr[1]
	
	# do some outlier trimming (1.5 IQR)
	tmpIQR<-IQR(dataVector)
	quant<-quantile(dataVector, c(0.25, 0.5, 0.75), type = 1)
	q1<-quant[1]
	q2<-quant[2]
	q3<-quant[3]
	tmpMedian<-median(dataVector)
	tmpTop<-tmpMedian+(1.5*tmpIQR)
	tmpBottom<-tmpMedian-(1.5*tmpIQR)
	
	dataVector<-subset(dataVector,dataVector>tmpBottom)
	dataVector<-subset(dataVector,dataVector<tmpTop)
	
	variance<-var(dataVector)
	stdev<-sd(dataVector)
	
	barData<-c(barData,variance)
	

	print(paste("iqr=",tmpIQR,"q1=",q1,"q2=",q2,"q3=",q3,"median=",tmpMedian,"iqrTop=",tmpTop,"iqrBottom=",tmpBottom,"variance=",variance,"stdev=",stdev))
	
	if(nchar(tmpName) > 50) {
		tmpName<-substr(tmpName,0,50)
	}
	nameData[i]<-paste("(",length(dataVector),") ",tmpName,sep="")
	nameData[i]<-paste(tmpName,sep="")
	
}

barData

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