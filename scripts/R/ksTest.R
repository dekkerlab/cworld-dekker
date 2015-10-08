options(bitmapType='cairo')

args <- commandArgs(TRUE)
wd<-args[1]
pvalueString<-args[2]
replicateMean<-as.numeric(args[3])
replicateStdev<-as.numeric(args[4])

pvalues<-unlist(strsplit(pvalueString,","))
pvalues<-as.numeric(pvalues)

ksTest<-ks.test(pvalues,pnorm,mean=replicateMean,sd=replicateStdev,alternative="two.sided")
ksTest$p.value
ksTest<-ks.test(pvalues,pnorm,mean=replicateMean,sd=replicateStdev,alternative="l")
ksTest$p.value
ksTest<-ks.test(pvalues,pnorm,mean=replicateMean,sd=replicateStdev,alternative="g")
ksTest$p.value


