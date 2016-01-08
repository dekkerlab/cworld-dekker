options(bitmapType='cairo')
library("ggplot2")
library("reshape")

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
	cat("\t",i,"\t",LETTERS[i],"\t",arg,"\n")
}

cat("\n")

allX <- vector()
myData <- list()
myLabel <- list()
nameData <- list()
for (i in 1:nFiles) {
	
	tmpFile<-fileArray[i]
	tmpName<-basename(tmpFile)
	tmpName<-unlist(strsplit(tmpName, ".", fixed = TRUE))[1]
	
	cat(paste("loading ",tmpName," ... ",sep=""))
	
	tmpData<-read.table(tmpFile,skip=skip,header=F,sep="\t",stringsAsFactors=FALSE,comment.char="#")
	dataVector=tmpData[,eval(dataColumn)]
	dataVector<-subset(dataVector,dataVector!="NaN")
	dataVector<-subset(dataVector,dataVector!="NA")
	dataVector<-subset(dataVector,dataVector!=".")
	
	dataVector<-as.numeric(dataVector)
	
	dataVector<-sort(dataVector)
	
	dat.df <- data.frame(signal = dataVector,label = rep(tmpName, length(dataVector)))
	
	allX<-c(allX,dataVector)
	myData[[i]]<-dat.df
	
	nameData[i]<-paste("(",length(dataVector),") ",tmpName,sep="")
	cat(paste("done","\n",sep=""))
}

cat("\n")

# sorting
cat(paste("sorting ... ",sep=""))
allX.length<-length(allX)
allX.sorted<-sort(allX)
cat(paste("done","\n",sep=""))

x.tmp.min<-floor(allX.sorted[ceiling(allX.length*0.001)])
x.tmp.max<-ceiling(allX.sorted[floor(allX.length*0.999)])

if(is.na(minX)) {
	minX<-x.tmp.min
}
if(is.na(maxX)) {
	maxX<-x.tmp.max
}

cat(paste("merging ... ",sep=""))
# merga all df
allData.df <- do.call("rbind", myData)
cat(paste("done","\n",sep=""))

cat("\n")

cat(paste("plotting","\n",sep=""))

# png file
png(file=paste(name,".png",sep=""),height=600,width=900)
ggplot(allData.df, aes(x = signal, fill = label)) + geom_density(alpha=0.25) + xlim(minX,maxX) + theme_bw() + theme(legend.position="bottom",legend.direction="vertical")+ geom_vline(xintercept = 0, colour="black", linetype = "longdash") +  ggtitle("Multi-Density Plot")
dev.off()

# pdf file
pdf(file=paste(name,".pdf",sep=''))
ggplot(allData.df, aes(x = signal, fill = label)) + geom_density(alpha=0.25) + xlim(minX,maxX) + theme_bw() + theme(legend.position="bottom",legend.direction="vertical")+ geom_vline(xintercept = 0, colour="black", linetype = "longdash") +  ggtitle("Multi-Density Plot")
dev.off()

#ggsave()
