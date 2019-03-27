rm(list=ls())
library(ggplot2)
library(preprocessCore)
library(matrixStats)

data=read.table("accelerometery_aggregate_phenotypes.continuous.txt",header=TRUE)
data[data==-1000]=NA

#UsualWalkingPace
data=subset(data,select=c("FID","IID","UsualWalkingPace"))

#binarize the data 
data[data<0]=NA
slow_vs_other=data 
slow_vs_other[slow_vs_other==1]=4
slow_vs_other[slow_vs_other==3]=1
slow_vs_other[slow_vs_other==2]=1
slow_vs_other[slow_vs_other==4]=2

brisk_vs_other=data 
brisk_vs_other[brisk_vs_other==2]=1
brisk_vs_other[brisk_vs_other==3]=2

slow_vs_other[slow_vs_other==NA]=-1000
brisk_vs_other[brisk_vs_other==NA]=-1000
write.table(slow_vs_other,file="slow_vs_other.UsualWalkingPace.txt",sep='\t',quote=FALSE,row.names=FALSE,col.names=TRUE)
write.table(brisk_vs_other,file="brisk_vs_other.UsualWalkingPace.txt",sep='\t',quote=FALSE,row.names=FALSE,col.names=TRUE)

