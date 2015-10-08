options(bitmapType='cairo')

args <- commandArgs(TRUE)
inputFile<-args[1]

jobName<-basename(inputFile)

tmpData<-read.table(inputFile,header=T,sep="\t")

pngfile<-paste(jobName,".png",sep='')
png(pngfile,height=800,width=800)

plot(tmpData$interactionDistance,tmpData$correlation,main="correlaton per distance",xlab="Genomic Distance (bp)",ylab="Correlation Coefficient",type="n",ylim=c(-1,1),yaxt="n")
abline(h=0,col="black",lty=2)
abline(v=4000000,col="black",lty=2)
axis(2,seq(from=-1,to=1,length.out=11))

lines(tmpData$interactionDistance,tmpData$correlation,col="black",lwd=3)

dev.off()