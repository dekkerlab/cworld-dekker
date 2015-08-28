args <- commandArgs(TRUE)
wd<-args[1]
valueString<-args[2]
replicateMean<-as.numeric(args[3])
replicateStdev<-as.numeric(args[4])

values<-unlist(strsplit(valueString,","))
values<-as.numeric(values)

subsetMean<-mean(values)
subsetSd<-sd(values)

standardError<-(replicateStdev/sqrt(length(values)))
zScore<-((subsetMean-replicateMean)/standardError)

pvalue<-2*(pnorm(-abs(subsetMean),mean=replicateMean,sd=replicateStdev))
pvalue<-2*(pnorm(-abs(zScore)))
pvalue
