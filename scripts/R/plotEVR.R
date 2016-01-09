options(bitmapType='cairo')

args <- commandArgs(TRUE)
cwd<-args[1]
inputFile<-args[2]
jobName<-args[3]

setwd(cwd)

data<-read.table(inputFile,header=T,sep="\t",comment.char="#")

dat.df <- data.frame(evr = data$evr,eigenvector = data$eigenvector)
# turn into percent
dat.df$evr<-dat.df$evr*100

pngfile<-paste(inputFile,".png",sep='')
png(pngfile,height=400,width=600)
    
plot(dat.df$eigenvector,dat.df$evr,main="Eigenvector explained variance ratio",xlab="eigenvector 1 .. n",ylab="explained variance",ylim=c(0,100))
lines(dat.df$eigenvector,dat.df$evr)

dev.off()
