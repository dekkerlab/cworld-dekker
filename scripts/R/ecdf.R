options(bitmapType='cairo')

ks.boot <- function(x, y, ...,
		alternative = c("two.sided", "less", "greater"),
		exact = NULL, nboots = 50) {
		
	# create progress bar
	pb <- txtProgressBar(min = 0, max = nboots, title="ks.boot", style = 3)
	
	alt <- match.arg(alternative)
	n <- length(x)
	D <- numeric(nboots)
	p <- numeric(nboots)
	 for(i in seq_len(nboots)){
		setTxtProgressBar(pb, i)
		idx <- sample(n, n, replace = TRUE)
		ks <- ks.test(x[idx], y, ..., alternative = alt, exact = exact)
		D[i] <- ks$statistic
		p[i] <- ks$p.value
	}
	list(D = mean(D), p.value = mean(p), nboots = nboots)
	
	close(pb)

}

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
	
	allX<-c(allX,dataVector)
	myData[[i]]<-dataVector
	
	nameData[i]<-paste("(",length(dataVector),") ",tmpName,sep="")
	cat(paste("done","\n",sep=""))
}

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

for (i1 in 1:nFiles) {
	y1<-myData[[i1]]
	name1<-nameData[[i1]]
	for (i2 in 1:nFiles) {
		
		if(i1 <= i2) {
			next
		}
		
		y2<-myData[[i2]]
		name2<-nameData[[i2]]
		
		#kstest<-ks.test(y1,y2,alternative="two.sided",exact=TRUE)
		#kstest.pval<-kstest$p.value
		
		#ksboot<-ks.boot(y1,y2,alternative="two.sided",nboots=100,exact=TRUE)
		#ksboot.pval<-ksboot$p.value
		
		#wilcox<-wilcox.test(y1,y2,paired=FALSE,alternative= "two.sided",exact=TRUE)
		#wilcox.pval<-wilcox$p.value
		
		#cat(paste("\tks.test\t",kstest.pval,"\n",sep=""))
		#cat(paste("\n",i1,"\t",name1,"\n",i2,"\t",name2,"\n",sep=""))
		#cat(paste("\tks.test\t",kstest.pval,"\twilcox\t",wilcox.pval,"\n",sep=""))
		
		#cat(paste("\tks.test\t",kstest.pval,"\tks.boot\t",ksboot.pval,"\twilcox\t",wilcox.pval,"\n",sep=""))
		
	}
}

cat("\n")

cat(paste("plotting","\n",sep=""))

# png file

pngfile<-paste(name,".png",sep='')
png(pngfile,height=800,width=800)

par(mar=c(5,5,5,5))

plot(0,main=paste("cumlative fraction"),ylab="cumlative fraction",xlab="x",xlim=c(minX,maxX),ylim=c(0,1),type="n",lwd=1,lty=1)
lineColors<-rainbow(nFiles,alpha=0.5)

for (i in 1:nFiles) {
	x<-myData[[i]]
	tmpName<-nameData[[i]]
	tmpColor=lineColors[i]
	
	lines(ecdf(x),col=tmpColor,lwd=1,lty=1,verticals="true")
}


legend("topleft", inset=0.01, legend=paste(nameData,sep=""), text.col=lineColors,cex=0.8)

dev.off()

# pdf file

pdffile<-paste(name,".pdf",sep='')
pdf(pdffile)

par(mar=c(5,5,5,5))

plot(0,main=paste("cumlative fraction"),ylab="cumlative fraction",xlab="x",xlim=c(minX,maxX),ylim=c(0,1),type="n",lwd=1,lty=1)

for (i in 1:nFiles) {
	x<-myData[[i]]
	tmpName<-nameData[[i]]
	tmpColor=lineColors[i]
	
	ecdf.d <- ecdf(x)
	
	## Construct a "linear interpolator" using approxfun:
	poly.ecdf <- with(environment(ecdf.d), approxfun(x,y))

	## 'test' it graphically, using curve():
	curve(poly.ecdf(x), add = TRUE, col=tmpColor,lwd=1,lty=1)
}

legend("topleft", inset=0.01, legend=paste(nameData,sep=""), text.col=lineColors,cex=0.8)

dev.off()

