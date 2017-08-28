rm(list=ls())
#check mutual information 

#post-quantile-normalization 
data=read.table("accelerometery_continuous_features_quantile_normalized_no_outliers.tsv",header=TRUE,sep='\t')
data[data==-1000]=NA
#attach(data)
library(infotheo)
library(ggplot2)
self_rep=c("DurationOfWalks","FrequencyStrenuousSportsLast4Weeks","FrequencyWalkingForPleasure",
           "NumberOfDaysModeratePhysicalActivity","DurationModerateActivity","NumberOfDaysVigorousPhysicalActivity",
           "DurationVigorousActivity","TimeSpentOutdoorsSummer","TimeSpentOutdoorsWinter","DurationWalkingForPleasure",
           "DurationStrenuousSports","NumberDaysWalked10Minutes","UsualWalkingPace")

measured=c("OverallAccelerationAverage","X6_7","X7_8","X12_13","X17_18","X18_19","StandardDeviationOfAcceleration","X9mg","X50mg") 

#perform binning! 
data$OverallAccelerationAverage=round(data$OverallAccelerationAverage/10)*10
data$X6_7=round(data$X6_7/10)*10
data$X7_8=round(data$X7_8/10)*10
data$X12_13=round(data$X12_13/10)*10
data$X17_18=round(data$X17_18/10)*10
data$X18_19=round(data$X18_19/10)*10
data$StandardDeviationOfAcceleration=round(data$StandardDeviationOfAcceleration/10)*10
data$X9mg=round(data$X9mg,0.1)
data$X50mg=round(data$X50mg,0.1)

data$DurationOfWalks=round(data$DurationOfWalks/10)*10
data$DurationModerateActivity=round(data$DurationModerateActivity/10)*10
data$DurationVigorousActivity=round(data$DurationVigorousActivity/10)*10
data$FrequencyStrenuousSportsLast4Weeks=round(data$FrequencyStrenuousSportsLast4Weeks/10)*10
data$FrequencyWalkingForPleasure=round(data$FrequencyWalkingForPleasure/10)*10

num_self_rep=length(self_rep)
num_measured=length(measured)
for(i in seq(1,num_measured)){
  f1=measured[i]
  for(j in seq(1,num_self_rep))
  {
    f2=self_rep[j]
    tmp=data.frame(data[[f1]],data[[f2]])
    names(tmp)=c(f1,f2)
    tmp=na.omit(tmp)
    mi=mutinformation(tmp[[f1]],tmp[[f2]])
    fname=paste(f1,'_vs_',f2,'.mutinfo.png',sep="")
    png(fname)
    smoothScatter(tmp[[f1]],
                  tmp[[f2]],
                  xlab=f1,
                  ylab=f2,
                  main=paste("Q-Normed ", "Mutual Info=",mi))
    dev.off() 
    
  }
}