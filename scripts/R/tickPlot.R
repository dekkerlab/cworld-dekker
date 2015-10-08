options(bitmapType='cairo')

args <- commandArgs(TRUE)
inputFile<-args[1]

inputFileName<-basename(inputFile)

myTable <- as.matrix(read.table(inputFile, header=TRUE, sep = "\t",row.names = 1,check.names=FALSE),check.names=FALSE)

pngfile<-paste(inputFileName,".png",sep='')
png(pngfile,height=800,width=800)

par(mar=c(10,5,5,10))

barplot(myTable,col=c("gray","orange","red"),xlab="",ylab="percent",las=2)

par(xpd=TRUE)
legend("topright",inset=c(-0.2,0),legend=c("nonelement:nonelement","element:nonelement","element:element"),fill=c("gray", "orange", "red"))

dev.off()
