options(bitmapType='cairo')

args <- commandArgs(TRUE)
inputFile<-args[1]

inputFileName<-basename(inputFile)

myTable <- read.table(inputFile, header=TRUE, sep = "\t",row.names = 1,check.names=FALSE)

for(group in colnames(myTable)) {
        tmpTable<-myTable[[group]]

        if( (max(tmpTable)==0) && (min(tmpTable)==0) ) {
                next
        }

        pngfile<-paste(inputFileName,".",group,".png",sep='')
        png(pngfile,height=800,width=1600)

        par(mar=c(10,5,5,10))

        tmpMatrix<-as.matrix(tmpTable,check.names=FALSE)
        dimnames(tmpMatrix) = list(rownames(myTable),c("subMatrix"))

        barplot(t(tmpMatrix),main=paste(group),ylab="sum of signal",xlab="",beside=TRUE,las=2)

        dev.off()

}

cis<-myTable$cis__cis+myTable$trans__cis
trans<-myTable$cis__trans+myTable$trans__trans
cis_trans_ratio<-cis/trans

pngfile<-paste(inputFileName,".cisTransRatio.png",sep='')
png(pngfile,height=800,width=1600)

par(mar=c(10,5,5,10))

tmpMatrix<-as.matrix(cis_trans_ratio,check.names=FALSE)
dimnames(tmpMatrix) = list(rownames(myTable))

barplot(t(tmpMatrix),main="cis/trans ratio",ylab="cis/trans",xlab="",beside=TRUE,las=2)

dev.off()