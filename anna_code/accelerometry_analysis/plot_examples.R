rm(list=ls())
library(wavelets)
ss<-function(x)
{
  sum(x^2)
}

#plotting some outliers 
outliers=c("2965229","2825432","5436759","1432711","2709188","5248604","5358474","5822979")
for(outlier in outliers){
  data=na.omit(as.numeric(as.matrix(read.table(paste("outliers/",outlier,"_90004_0_0.csv",sep="")))))
  tran=dwt(data,fast=TRUE,n.level=13)
  png(paste(outlier,".png",sep=""))
  plot.dwt(tran)
  dev.off() 
}

#raw examples that are not outliers 
examples=c("1880013","1843818","1650596","1563156","1352287","1313911","1395004","1886476","1129255","1496521")
for(example in examples){
  data=na.omit(as.numeric(as.matrix(read.table(paste("examples/",example,"_90004_0_0.csv",sep="")))))
  tran=dwt(data,fast=TRUE,n.level=13)
  png(paste(example,".png",sep=""))
  plot.dwt(tran)
  dev.off() 
}
