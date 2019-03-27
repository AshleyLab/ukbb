rm(list=ls())
library(ggplot2)
library(matrixStats)
library(preprocessCore)
library(reshape2)
data=read.table("../accelerometry_analysis/averages/subject_acceleration_interval_averages.txt",header=TRUE,row.names=1)
#use a cutoff of 200 mg's 
#use a cutoff of 15 mg's (based on visual analysis of the skew in the histogram)
too_high1=which(data[,1]>10)
data[too_high1,1]=NA
too_high5=which(data[,5]>25)
data[too_high5,5]=NA
too_high8=which(data[,8]>7.5)
data[too_high8,8]=NA
for(col in c(2,3,4,5,6,7,9))
{
  too_high=which(data[,col]>200)
  data[too_high,col]=NA
}
melted=melt(data)
melted$group[melted$variable=="X21_29"]=1
melted$group[melted$variable=="X29_45"]=1
melted$group[melted$variable=="X33_41"]=2
melted$group[melted$variable=="X41_48"]=2
melted$group[melted$variable=="X48_57"]=2
melted$group[melted$variable=="X12_18"]=3
melted$group[melted$variable=="X18_24"]=3
melted$group[melted$variable=="X24_30"]=3
melted$group[melted$variable=="X30_36"]=3
melted$group=factor(melted$group)
ggplot(melted,aes(melted$variable,melted$value,fill=melted$group))+
  geom_boxplot()+
  xlab("Time Intervals")+
  ylab("Subject mean acceleration (mg)")+
  scale_x_discrete(breaks=c("X21_29", "X29_45","X33_41","X41_48","X48_57","X12_18","X18_24","X24_30","X30_36"),
                   labels=c("9pm-5am (night)","5am - 9pm (day)","9am-5pm (work)","5pm-12am (evening)","12am-9am (night + morn.)","12pm-6pm (afternoon)","6pm-12am (evening)","12am-6am (night)","6am-12pm (morning)"))+
  theme_bw(20)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
data$FID=as.numeric(rownames(data))
#get residuals
covar=data.frame(read.table('covariates.txt',header=TRUE,sep='\t'))
covar[covar==-1000]=NA
covar=covar[order(covar$FID),]
covar=covar[4:nrow(covar),]
data=data[order(data$FID),]
covar$Sex=factor(covar$Sex)
covar$f.batch=factor(covar$f.batch)
residuals_continuous=matrix(nrow=nrow(data),ncol=ncol(data))

merged=merge(data,covar,by="FID")
residuals_continuous[,1]=merged[,1]
for(col in seq(2,10))
{
merged$Y=merged[,col]
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
                                                +f.batch,data=merged,na.action=na.exclude),na.action=na.exclude))
residuals_continuous[,col]=residuals
}
residuals_continuous=as.data.frame(residuals_continuous)
names(residuals_continuous)=names(merged)[1:10]

#quantile normalize the residuals 
for(col in seq(3,ncol(residuals_continuous)))
{
residuals_continuous[,col]=normalize.quantiles(as.matrix(residuals_continuous[,col]))
}
residuals_continuous$IID=residuals_continuous$FID 
residuals_continuous[is.na(residuals_continuous)]=-1000
residuals_continuous=residuals_continuous[c(1,11,2,3,4,5,6,7,8,9,10)]
names(residuals_continuous)=c("FID","IID","X9pm_5am","X5am_9pm","X9am_5pm","X5pm_12am","X12am_9am","X12pm_6pm","X6pm_12am","X12am_6am","X6am_12pm")
write.table(residuals_continuous,file="acceleration.average.custom.intervals.txt",sep='\t',quote=FALSE,row.names=FALSE,col.names=TRUE)
