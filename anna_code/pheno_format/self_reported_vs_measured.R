rm(list=ls())
#check correlation between self-reported and measured features 

#post-quantile-normalization 
data=read.table("accelerometery_continuous_features_quantile_normalized_no_outliers.tsv",header=TRUE,sep='\t')
#pre-quantile normalization 
#data=read.table("accelerometery_aggregate_phenotypes.continuous.txt",header=TRUE,sep='\t')
data[data==-1000]=NA
attach(data)
library(ggplot2)
self_rep=c("DurationOfWalks","FrequencyStrenuousSportsLast4Weeks","FrequencyWalkingForPleasure",
           "NumberOfDaysModeratePhysicalActivity","DurationModerateActivity","NumberOfDaysVigorousPhysicalActivity",
           "DurationVigorousActivity","TimeSpentOutdoorsSummer","TimeSpentOutdoorsWinter","DurationWalkingForPleasure",
           "DurationStrenuousSports","NumberDaysWalked10Minutes","UsualWalkingPace")

measured=c("OverallAccelerationAverage","X6_7","X7_8","X12_13","X17_18","X18_19","StandardDeviationOfAcceleration","X9mg","X50mg") 
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
    corval=cor(tmp[[f1]],tmp[[f2]],method="spearman")
    fname=paste(f1,'_vs_',f2,'.png',sep="")
    png(fname)
    smoothScatter(tmp[[f1]],
                  tmp[[f2]],
                  xlab=f1,
                  ylab=f2,
                  main=paste("Q-Normed ", "spearman=",corval))
   dev.off() 
    
  }
}