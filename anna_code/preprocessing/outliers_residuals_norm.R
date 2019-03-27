rm(list=ls())
library(ggplot2)
library(preprocessCore)
library(matrixStats)

data=read.table("accelerometery_aggregate_phenotypes.continuous.txt",header=TRUE)
data[data==-1000]=NA
#remove any outliers more than 3 st dev away from mean 
d_mean=colMeans(data,na.rm=TRUE)
d_sd=colSds(as.matrix(data),na.rm=TRUE )
upper_bound=d_mean+3*d_sd
lower_bound=d_mean-3*d_sd 
for(col in seq(3,ncol(data)))
{
  cur_upper_bound=upper_bound[col]
  cur_lower_bound=lower_bound[col]
  to_truncate_upper=which(data[,col]>cur_upper_bound)
  data[to_truncate_upper,col]=cur_upper_bound 
  to_truncate_lower=which(data[,col]<cur_lower_bound) 
  data[to_truncate_lower,col]=cur_lower_bound
}
data[data==NA]=-1000
write.table(data,file="accelerometry_aggregate_phenotypes.continuous.no_outliers.txt",sep='\t',quote=FALSE,row.names=FALSE,col.names=TRUE)
data[data==-1000]=NA
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
write.table(residuals_continuous,file="accelerometry_aggregate_phenotypes.continuous.no_outliers.residuals.qnorm.txt",sep='\t',quote=FALSE,row.names=FALSE,col.names=TRUE)
