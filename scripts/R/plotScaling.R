minor.ticks.axis <- function(ax,n,t.ratio=0.5,mn,mx){

  lims <- par("usr")
  if(ax %in%c(1,3)) lims <- lims[1:2] else lims[3:4]

  major.ticks <- seq(mn,mx,1)

  labels <- sapply(major.ticks,function(i)
            as.expression(bquote(10^ .(i)))
          )
  axis(ax,at=major.ticks,labels=labels)

  n <- n+2
  minors <- log10(pretty(10^major.ticks[1:2],n))-major.ticks[1]
  minors <- minors[-c(1,n)]

  minor.ticks = c(outer(minors,major.ticks,`+`))
  minor.ticks <- minor.ticks[minor.ticks > mn & minor.ticks < mx]

  axis(ax,at=minor.ticks,tcl=par("tcl")*t.ratio,labels=FALSE)
}

args <- commandArgs(TRUE)
dir<-args[1]
inputFiles<-args[2]
inputNames<-args[3]
averageTransSignals<-args[4]
jobName<-args[5]

setwd(dir)

fileArray<-unlist(strsplit(inputFiles,","))
nameArray<-unlist(strsplit(inputNames,","))
averageTransArray<-unlist(strsplit(averageTransSignals,","))
nFiles<-length(fileArray)

allX<-list()
allY<-list()
allData<-data.frame()

for (i in 1:nFiles) {
	tmpFile<-fileArray[i]

	tmpLoessData<-read.table(tmpFile,header=T,sep="\t",stringsAsFactors=FALSE)
	
	# remove 0s
	tmpLoessData<-subset(tmpLoessData,tmpLoessData$realInteractionDistance>0)
	tmpLoessData<-subset(tmpLoessData,tmpLoessData$loessExpectedValue>0)
	tmpLoessData<-subset(tmpLoessData,tmpLoessData$loessExpectedStdev>0)
	tmpLoessData<-subset(tmpLoessData,tmpLoessData$observedSignal>0)
	
	# remove NAs
	tmpLoessData<-subset(tmpLoessData,tmpLoessData$realInteractionDistance!="NA")
	tmpLoessData<-subset(tmpLoessData,tmpLoessData$loessExpectedValue!="NA")
	tmpLoessData<-subset(tmpLoessData,tmpLoessData$loessExpectedStdev!="NA")
	tmpLoessData<-subset(tmpLoessData,tmpLoessData$observedSignal!="NA")

	# log transform
	tmpLoessData$realInteractionDistance<-log10(tmpLoessData$realInteractionDistance)
	tmpLoessData$loessExpectedValue<-log10(tmpLoessData$loessExpectedValue)
	tmpLoessData$loessExpectedStdev<-log10(tmpLoessData$loessExpectedStdev)
	tmpLoessData$observedSignal<-log10(tmpLoessData$observedSignal)
	
	# subset all interactions > 10^3 (1000bp) or larger only
	tmpLoessData<-subset(tmpLoessData,tmpLoessData$realInteractionDistance>3)
	
	nrow(tmpLoessData)
	allData<-rbind(allData,tmpLoessData)
	
}

allX<-c(allData$realInteractionDistance,allData$realInteractionDistance)
allY<-c(allData$loessExpectedValue,allData$observedSignal)

ymin<-floor(min(allY))
ymax<-ceiling(max(allY))

xmin<-min(allX)
xmax<-max(allX)

#500kb - 7MB
minDistance<-log10(500000)
maxDistance<-log10(7000000)

distanceSpan<-(maxDistance-minDistance)
if(minDistance < xmin) {
	xmin<-(minDistance-(distanceSpan*0.25))
}
if(maxDistance > xmax) {
	xmax<-(maxDistance+(distanceSpan*0.25))
}

xmin<-floor(xmin)
xmax<-ceiling(xmax)

# now plot scaling data
pngfile<-paste(jobName,".png",sep='')
png(pngfile,height=800,width=800)

# plot all the data - setup area
plot(allX,allY,xlab="genomic distance",ylab="interaction score",main=paste("CData Scaling Plot","\n",jobName,sep=""),type="n",ylim=c(ymin,ymax),xlim=c(xmin,xmax),xaxt="n",yaxt="n")

if(nFiles > 1) {
	pointColors<-rainbow(nFiles,alpha=0.05)
	lineColors<-rainbow(nFiles,alpha=0.5)
} else {
	pointColors<-rgb(0.5,0.5,0.5,0.05)
	lineColors<-rgb(0,0,0,0.8)
}

slopeArray<-list()

for (i in 1:nFiles) {
	tmpFile<-fileArray[i]
	
	averageTransSignal<-as.numeric(averageTransArray[i])
	
	tmpLoessData<-read.table(tmpFile,header=T,sep="\t",stringsAsFactors=FALSE)	
	
	# remove 0s
	tmpLoessData<-subset(tmpLoessData,tmpLoessData$realInteractionDistance>0)
	tmpLoessData<-subset(tmpLoessData,tmpLoessData$loessExpectedValue>0)
	tmpLoessData<-subset(tmpLoessData,tmpLoessData$loessExpectedStdev>0)
	tmpLoessData<-subset(tmpLoessData,tmpLoessData$observedSignal>0)
	
	# remove NAs
	tmpLoessData<-subset(tmpLoessData,tmpLoessData$realInteractionDistance!="NA")
	tmpLoessData<-subset(tmpLoessData,tmpLoessData$loessExpectedValue!="NA")
	tmpLoessData<-subset(tmpLoessData,tmpLoessData$loessExpectedStdev!="NA")
	tmpLoessData<-subset(tmpLoessData,tmpLoessData$observedSignal!="NA")
	
	nDataPoints<-nrow(tmpLoessData)
	if(nDataPoints == 0) {
		next
	}
	
	# log transform 
	tmpLoessData$realInteractionDistance<-log10(tmpLoessData$realInteractionDistance)
	tmpLoessData$loessExpectedValue<-log10(tmpLoessData$loessExpectedValue)
	tmpLoessData$loessExpectedStdev<-log10(tmpLoessData$loessExpectedStdev)
	tmpLoessData$observedSignal<-log10(tmpLoessData$observedSignal)
	
	# subset all interactions > 10^3 (1000bp) or larger only
	tmpLoessData<-subset(tmpLoessData,tmpLoessData$realInteractionDistance>3)
	
	tmpLoessData<-tmpLoessData[with(tmpLoessData, order(realInteractionDistance)), ]	

	# draw the scaling line
	points(tmpLoessData$realInteractionDistance,tmpLoessData$observedSignal,col=pointColors[i],cex=0.5)	
	lines(tmpLoessData$realInteractionDistance,tmpLoessData$loessExpectedValue,col=lineColors[i],lwd=4)	
	
	# draw the average trans signal
	if(!is.na(averageTransSignal)) {
		abline(h=log10(averageTransSignal),col=lineColors[i],lwd=3,lty=2)
	}
	
	tmpLoessData<-subset(tmpLoessData,tmpLoessData$realInteractionDistance>=minDistance & tmpLoessData$realInteractionDistance<=maxDistance)
	
	y <- tmpLoessData$observedSignal
	x <- tmpLoessData$realInteractionDistance
	f <- lm(y~x)
	X<-c(minDistance-(distanceSpan*0.5),maxDistance+(distanceSpan*0.5))
	Y<-predict(f, newdata=data.frame(x=X))
	lines(x=X, y=Y,col=lineColors[i],lwd=2,lty=2)
	slope<-f$coefficients[[2]]
	slope<-round(slope, digits = 4)
	
	slopeArray[i]<-slope
}

for(i in 1:nFiles) {
	name<-nameArray[i]
	slope<-slopeArray[i]
	#if(nchar(name) > 40) {
	#	name<-substr(name,0,40)
	#}
	label<-paste(name," (",slope,")",sep="")
	
	nameArray[i]<-label
}

rect(minDistance, ymin, maxDistance, ymax, col=rgb(0,1,0,0.1),border=rgb(0,1,0,0.1))

minor.ticks.axis(1,9,mn=xmin,mx=xmax)
minor.ticks.axis(2,9,mn=ymin,mx=ymax)

legend("topright", inset=0.01, legend=paste(nameArray,sep=""), text.col=lineColors,cex=0.65)

dev.off()





# now plot scaling data
pdffile<-paste(jobName,".pdf",sep='')
pdf(pdffile)

# plot all the data - setup area
plot(allX,allY,xlab="genomic distance",ylab="interaction score",main=paste("CData Scaling Plot","\n",jobName,sep=""),type="n",ylim=c(ymin,ymax),xlim=c(xmin,xmax),xaxt="n",yaxt="n")

if(nFiles > 1) {
	pointColors<-rainbow(nFiles,alpha=0.05)
	lineColors<-rainbow(nFiles,alpha=0.8)
} else {
	pointColors<-rgb(0.5,0.5,0.5,0.05)
	lineColors<-rgb(0,0,0,0.8)
}

slopeArray<-list()

for (i in 1:nFiles) {
	tmpFile<-fileArray[i]
	
	averageTransSignal<-as.numeric(averageTransArray[i])
	
	tmpLoessData<-read.table(tmpFile,header=T,sep="\t",stringsAsFactors=FALSE)	
	
	# remove 0s
	tmpLoessData<-subset(tmpLoessData,tmpLoessData$realInteractionDistance>0)
	tmpLoessData<-subset(tmpLoessData,tmpLoessData$loessExpectedValue>0)
	tmpLoessData<-subset(tmpLoessData,tmpLoessData$loessExpectedStdev>0)
	tmpLoessData<-subset(tmpLoessData,tmpLoessData$observedSignal>0)
	
	# remove NAs
	tmpLoessData<-subset(tmpLoessData,tmpLoessData$realInteractionDistance!="NA")
	tmpLoessData<-subset(tmpLoessData,tmpLoessData$loessExpectedValue!="NA")
	tmpLoessData<-subset(tmpLoessData,tmpLoessData$loessExpectedStdev!="NA")
	tmpLoessData<-subset(tmpLoessData,tmpLoessData$observedSignal!="NA")
	
	nDataPoints<-nrow(tmpLoessData)
	if(nDataPoints == 0) {
		next
	}
	
	# log transform 
	tmpLoessData$realInteractionDistance<-log10(tmpLoessData$realInteractionDistance)
	tmpLoessData$loessExpectedValue<-log10(tmpLoessData$loessExpectedValue)
	tmpLoessData$loessExpectedStdev<-log10(tmpLoessData$loessExpectedStdev)
	tmpLoessData$observedSignal<-log10(tmpLoessData$observedSignal)
	
	# subset all interactions > 10^3 (1000bp) or larger only
	tmpLoessData<-subset(tmpLoessData,tmpLoessData$realInteractionDistance>3)
	
	tmpLoessData<-tmpLoessData[with(tmpLoessData, order(realInteractionDistance)), ]	

	# draw the scaling line
	points(tmpLoessData$realInteractionDistance,tmpLoessData$observedSignal,col=pointColors[i],cex=0.5)	
	lines(tmpLoessData$realInteractionDistance,tmpLoessData$loessExpectedValue,col=lineColors[i],lwd=4)	
	
	# draw the average trans signal
	if(!is.na(averageTransSignal)) {
		abline(h=log10(averageTransSignal),col=lineColors[i],lwd=3,lty=2)
	}
	
	tmpLoessData<-subset(tmpLoessData,tmpLoessData$realInteractionDistance>=minDistance & tmpLoessData$realInteractionDistance<=maxDistance)
	
	y <- tmpLoessData$observedSignal
	x <- tmpLoessData$realInteractionDistance
	f <- lm(y~x)
	X<-c(minDistance-(distanceSpan*0.5),maxDistance+(distanceSpan*0.5))
	Y<-predict(f, newdata=data.frame(x=X))
	lines(x=X, y=Y,col=lineColors[i],lwd=2,lty=2)
	slope<-f$coefficients[[2]]
	slope<-round(slope, digits = 4)
	
	slopeArray[i]<-slope
}

for(i in 1:nFiles) {
	name<-nameArray[i]
	slope<-slopeArray[i]
	if(nchar(name) > 40) {
		name<-substr(name,0,40)
	}
	label<-paste(name," (",slope,")",sep="")
	
	nameArray[i]<-label
}

rect(minDistance, ymin, maxDistance, ymax, col=rgb(0,1,0,0.1),border=rgb(0,1,0,0.1))

minor.ticks.axis(1,9,mn=xmin,mx=xmax)
minor.ticks.axis(2,9,mn=ymin,mx=ymax)

legend("topright", inset=0.01, legend=paste(nameArray,sep=""), text.col=lineColors,cex=0.8)
