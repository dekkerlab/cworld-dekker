# http://www.burns-stat.com/pages/Working/perfmeasrandport.pdf
args <- commandArgs(TRUE)
wd<-args[1]
pvalueString<-args[2]

pvalues<-unlist(strsplit(pvalueString,","))
pvalues<-as.numeric(pvalues)

stouffers.pvalue<-pnorm(sum(qnorm(pvalues)) / sqrt(length(pvalues)))

print(paste("stouffers=",stouffers.pvalue,sep=""))


