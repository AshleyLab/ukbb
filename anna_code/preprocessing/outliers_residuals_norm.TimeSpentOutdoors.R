rm(list=ls())
library(ggplot2)
library(preprocessCore)
library(matrixStats)

data=read.table("accelerometery_aggregate_phenotypes.continuous.txt",header=TRUE)
data[data==-1000]=NA

#TimeSpentOutdoors 
timeOutdoors=subset(data,select=c("FID","IID","TimeSpentOutdoorsSummer","TimeSpentOutdoorsWinter"))
ggplot(timeOutdoors,aes(x=TimeSpentOutdoorsSummer))+
  geom_histogram(stat="bin")

#set values of -10 to 0.5 (indicates less than 1 hour per day)
timeOutdoors[timeOutdoors==-10]=0.5 
ggplot(timeOutdoors,aes(x=TimeSpentOutdoorsSummer))+
  geom_histogram('binwdith'=1)

#remove values of -3 (Prefer not to answer) and -1 (Do not know)
timeOutdoors[timeOutdoors==-3]=NA
timeOutdoors[timeOutdoors==-1]=NA

ggplot(timeOutdoors,aes(x=TimeSpentOutdoorsSummer))+
  geom_histogram('binwdith'=2)

#remove any outliers more than 3 st dev away from mean 
d_mean=colMeans(timeOutdoors,na.rm=TRUE)
d_sd=colSds(as.matrix(timeOutdoors),na.rm=TRUE )
upper_bound=d_mean+3*d_sd
lower_bound=d_mean-3*d_sd 
for(col in seq(3,ncol(timeOutdoors)))
{
  cur_upper_bound=upper_bound[col]
  cur_lower_bound=lower_bound[col]
  to_truncate_upper=which(timeOutdoors[,col]>cur_upper_bound)
  timeOutdoors[to_truncate_upper,col]=cur_upper_bound 
  to_truncate_lower=which(timeOutdoors[,col]<cur_lower_bound) 
  timeOutdoors[to_truncate_lower,col]=cur_lower_bound
}
#get residuals
covar=data.frame(read.table('covariates.txt',header=TRUE,sep='\t'))
covar[covar==-1000]=NA
covar=covar[order(covar$FID),]
covar=covar[4:nrow(covar),]
timeOutdoors=timeOutdoors[order(timeOutdoors$FID),]
covar$Sex=factor(covar$Sex)
covar$f.batch=factor(covar$f.batch)
residuals_continuous=matrix(nrow=nrow(timeOutdoors),ncol=ncol(timeOutdoors))
residuals_continuous[,1]=timeOutdoors[,1]
residuals_continuous[,2]=timeOutdoors[,2]
for(col in seq(3,ncol(timeOutdoors)))
{
covar$Y=timeOutdoors[,col]
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
names(residuals_continuous)=names(timeOutdoors)

#quantile normalize the residuals 
for(col in seq(3,ncol(residuals_continuous)))
{
residuals_continuous[,col]=normalize.quantiles(as.matrix(residuals_continuous[,col]))
}
residuals_continuous[residuals_continuous==NA]=-1000
write.table(residuals_continuous,file="TimeOutdoors.no_outliers.residuals.qnorm.txt",sep='\t',quote=FALSE,row.names=FALSE,col.names=TRUE)
