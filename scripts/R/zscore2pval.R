options(bitmapType='cairo')

library(fBasics)
library(qvalue)

args <- commandArgs(TRUE)
inputFile<-args[1]
inputName<-basename(inputFile)

loessData <- read.table(inputFile,header=TRUE)

print('attemping to model null ZSCOREs')
print('Assuming all negative zScores are true NULLs')
quant <- quantile(loessData[,"zScore"][loessData[,"zScore"]<0], c(0,.01, .02, .05),na.rm=TRUE)
print(paste('1%-ile of negative zScores is ',quant[2]))
thresh <- quantile(loessData[,"zScore"][loessData[,"zScore"]<0], .01,na.rm=TRUE)
print(paste('Using ',thresh,'<= ZSCOREs <=',-thresh,' to model true NULLs'))

threshz <- loessData[,"zScore"][loessData[,"zScore"] <= -thresh]
shiftz <- threshz[threshz > thresh] - thresh
wfit <- fitdistr(shiftz, "weibull")
shape <- wfit$estimate["shape"]
scale <- wfit$estimate["scale"]

wfit

par(mfrow=c(1,2))

png(paste(inputName,'.fit.all.png',sep=''), height=400, width=600)
hist(threshz, breaks = 1000, probability = T, main = "z-scores | Weibull fit to null",xlab = 'z-score')
lines(seq(0,4,length=100)+thresh, dweibull(seq(0,4,length=100), shape=shape, scale=scale), col = 2, lwd = 2)
dev.off()

png(paste(inputName,'.fit.lt10.png',sep=''), height=400, width=600)
hist(loessData[,"zScore"][loessData[,"zScore"]<10], breaks = 1000, probability = T, main = "z-scores < 10 |  Weibull fit to null",xlab = 'z-score')
lines(seq(0,10,length=500)+thresh, dweibull(seq(0,10,length=500), shape=shape, scale=scale), col = 2, lwd = 2)
abline(v=-thresh, col = 4, lwd = 2)
dev.off()

wq <- qweibull(ppoints(length(shiftz)), shape=shape, scale=scale)
wqlin <- lm(z~wq, data.frame(z=sort(shiftz), wq=wq))
summary(wqlin)

wp <- rep(1, nrow(loessData))
wp[loessData[,"zScore"] > thresh] <- 1-pweibull(loessData[,"zScore"][loessData[,"zScore"] > thresh] - thresh, shape=shape, scale=scale)
qvals <- qvalue(wp, pi0.method = 'bootstrap')
write.table(data.frame(loessData$yHeaderName,loessData$xHeaderName,loessData$interactionDistance,loessData$observedSignal,loessData$loessExpectedValue,loessData$loessExpectedStdev,loessData$zScore,wp, qvals$q), file = paste(inputName,'.zscore-pval-qval.txt',sep=''), quote = F,sep = '\t', col.names = F, row.names = F)
