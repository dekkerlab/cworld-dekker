options(bitmapType='cairo')

args <- commandArgs(TRUE)
dir<-args[1]
inputFile<-args[2]
span<-as.numeric(args[3])

inputFileArray<-unlist(strsplit(inputFile, "\\/"))
inputFileName<-inputFileArray[length(inputFileArray)]

wd<-paste(dir,sep='')
setwd(wd)

tmpData<-read.table(inputFile,header=F,sep="\t")

tmpData<-subset(tmpData,tmpData[,1]!="NA")
tmpData<-subset(tmpData,tmpData[,2]!="NA")

label<-tmpData[,1]
x<-tmpData[,2]
y<-tmpData[,3]

tmpData.size=nrow(tmpData)

loessSpan<-tmpData.size*span
while(loessSpan < 10) { 
	span<-span*2
	print(paste("increasing loess span -> ",span))
    loessSpan<-tmpData.size*span    
}

loessData <- loess(y~x,span=span,family="symmetric")
smoothedData <- predict(loessData)

write.table(data.frame(tmpData[,1],tmpData[,2],tmpData[,3],smoothedData), file = paste(inputFileName,'.smoothed',sep=''), quote = F,sep = '\t', col.names = F, row.names = F)

