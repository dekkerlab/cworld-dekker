options(bitmapType='cairo')

args <- commandArgs(TRUE)
inputFile<-args[1]
correlationMode<-args[2]

if(!exists(correlationMode)) {
	correlationMode<-"pearson"
}

if((file.info(inputFile)[[1]][1]) > 0) {
	
	tmpData<-read.table(inputFile,header=TRUE,sep="\t")
	
	r<-cor(tmpData$score_1,tmpData$score_2,method=correlationMode);

	cat(r)
	
} else {
	cat("NA")
}