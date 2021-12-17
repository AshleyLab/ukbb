rm(list=ls())
#calculate correlation between exposure and outcome values. 
load("survival_data.RData")
exposures=read.table("mr_exposures.euro.txt",header=TRUE,sep='\t')
exposures=exposures[order(exposures$Subject),]
outcomes=read.table("mr_outcome.euro.txt",header=TRUE,sep='\t')
outcomes=outcomes[order(outcomes$Subject),]
subjects=outcomes$Subject
status=status[subjects]
time_to_death_years=time_to_death_years[subjects]
cvd=outcomes$CvdStatus
out_df=data.frame(time_to_death_years,cvd)
rownames(out_df)=subjects 
rownames(exposures)=subjects
exposures$Subject=NULL
for(exposure in names(exposures)){
  for(outcomes in names(out_df)){
    browser() 
  }
}