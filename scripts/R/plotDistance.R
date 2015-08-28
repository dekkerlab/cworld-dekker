minor.ticks.axis <- function(ax,n,t.ratio=0.5,mn,mx){

  lims <- par("usr")
  if(ax %in%c(1,3)) lims <- lims[1:2] else lims[3:4]

  major.ticks <- seq(floor(mn),ceiling(mx),1)

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
inputFile<-args[2]
jobName<-args[3]
logTransform<-as.numeric(args[4])

wd<-paste(dir,sep='')
setwd(wd)

tmpData<-read.table(inputFile,header=T,sep="\t")

tmpData<-subset(tmpData,tmpData$interactionDistance!="NA")

if(logTransform == 10) {
	# remove 0/NA
	tmpData<-subset(tmpData,tmpData$interactionDistance>0)
	# log transform
	tmpData$interactionDistance<-log10(tmpData$interactionDistance)
}

pngfile<-paste(jobName,".png",sep='')
png(pngfile,height=800,width=800)

xmin<-min(tmpData$interactionDistance)
xmax<-max(tmpData$interactionDistance)
xmin<-0

plot(tmpData$interactionDistance,tmpData$cumulativePercent,ylim=c(0,100),main=paste(jobName,"\n","cumulative reads/distance",sep=""),xlab="Genomic Distance (bp)",ylab="Cumulative Read Percent",type="n",xaxt="n",yaxt="n")
lines(tmpData$interactionDistance,tmpData$cumulativePercent,col="black",lwd=3)

axis(2,seq(from=0,to=100,by=10))

if(logTransform == 10) {
	minor.ticks.axis(1,9,mn=xmin,mx=xmax)
} else {
	axis(1,seq(from=xmin,to=xmax,length.out=5))	
}

dev.off()