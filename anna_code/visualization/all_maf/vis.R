rm(list=ls())
library(qqman)
args <- commandArgs(TRUE)
fname=args[1]
data=read.table(fname,header=TRUE,sep='\t')
data$P[data$P==0]=1e-100
manh_name=paste(fname,'.manh.png','_')
qq_name=paste(fname,'.qq.png','_')
png(manh_name)
manhattan(data,main=fname)
dev.off()
png(qq_name)
qq(data$P,main=fname)
dev.off()
