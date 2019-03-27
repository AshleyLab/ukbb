rm(list=ls())
library(ggplot2)
library(preprocessCore)
library(matrixStats)

data=read.table("accelerometery_aggregate_phenotypes.continuous.txt",header=TRUE)
data[data==-1000]=NA

#hourly acceleration averages 
data=subset(data,select=c("FID","IID","X0_1","X1_2","X2_3","X3_4","X4_5","X5_6"))

#use a cutoff of 15 mg's (based on visual analysis of the skew in the histogram)
for(col in seq(3,ncol(data)))
{
  too_high=which(data[,col]>15)
  data[too_high,col]=NA
}
#get residuals
covar=data.frame(read.table('covariates.txt',header=TRUE,sep='\t'))
covar[covar==-1000]=NA
covar=covar[order(covar$FID),]
covar=covar[4:nrow(covar),]
data=data[order(data$FID),]
covar$Sex=factor(covar$Sex)
covar$f.batch=factor(covar$f.batch)
residuals_continuous=matrix(nrow=nrow(data),ncol=ncol(data))
residuals_continuous[,1]=data[,1]
residuals_continuous[,2]=data[,2]
for(col in seq(3,ncol(data)))
{
covar$Y=data[,col]
residuals=as.vector(residuals(lm(Y ~ Sex
                                                +YearOfBirth
                                                +PC1
                                                +PC2
                                                +PC3
                                                +PC4
                                                +PC5
                                                +PC6
                                                +PC7
                                                +PC8
                                                +PC9
                                                +PC10
                                                +f.batch,data=covar,na.action=na.exclude),na.action=na.exclude))
residuals_continuous[,col]=residuals
}
residuals_continuous=as.data.frame(residuals_continuous)
names(residuals_continuous)=names(data)

#quantile normalize the residuals 
for(col in seq(3,ncol(residuals_continuous)))
{
residuals_continuous[,col]=normalize.quantiles(as.matrix(residuals_continuous[,col]))
}
residuals_continuous[residuals_continuous==NA]=-1000
write.table(residuals_continuous,file="hourly.acceleration.night.txt",sep='\t',quote=FALSE,row.names=FALSE,col.names=TRUE)
