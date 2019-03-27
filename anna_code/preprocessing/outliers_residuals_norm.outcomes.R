rm(list=ls())
library(ggplot2)
library(preprocessCore)
library(matrixStats)

data=read.table("mr_outcome.txt",header=TRUE)
data=data[order(data$Subject),]
categorical_data=subset(data,select=c("Subject","DeathStatus","CvdStatus","AliveAt70","AliveAt65"))
categorical_data[categorical_data=="NA"]=-1000 
categorical_data[categorical_data==1]=2 
categorical_data[categorical_data==0]=1 
write.table(categorical_data,file="mr_outcome.categorical.txt",sep='\t',quote=FALSE,row.names=FALSE,col.names=TRUE)

data=subset(data,select=c("Subject","DeathAge","AssessedAge"))
#remove any outliers more than 3 st dev away from mean 
d_mean=colMeans(data,na.rm=TRUE)
d_sd=colSds(as.matrix(data),na.rm=TRUE )
upper_bound=d_mean+3*d_sd
lower_bound=d_mean-3*d_sd 
for(col in seq(2,ncol(data)))
{
  cur_upper_bound=upper_bound[col]
  cur_lower_bound=lower_bound[col]
  to_truncate_upper=which(data[,col]>cur_upper_bound)
  data[to_truncate_upper,col]=cur_upper_bound 
  to_truncate_lower=which(data[,col]<cur_lower_bound) 
  data[to_truncate_lower,col]=cur_lower_bound
}
covar=data.frame(read.table('covariates.txt',header=TRUE,sep='\t'))
covar[covar==-1000]=NA
covar=covar[order(covar$FID),]
covar=covar[4:nrow(covar),]
data=data[order(data$Subject),]
covar$Sex=factor(covar$Sex)
covar$f.batch=factor(covar$f.batch)
residuals_continuous=matrix(nrow=nrow(data),ncol=ncol(data))
residuals_continuous[,1]=data[,1]
for(col in seq(2,ncol(data)))
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
for(col in seq(2,ncol(residuals_continuous)))
{
residuals_continuous[,col]=normalize.quantiles(as.matrix(residuals_continuous[,col]))
}
residuals_continuous[residuals_continuous==NA]=-1000
write.table(residuals_continuous,file="mr_outcome.continuous.no_outliers.residuals.qnorm.txt",sep='\t',quote=FALSE,row.names=FALSE,col.names=TRUE)
