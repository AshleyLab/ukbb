rm(list=ls())
library('MendelianRandomization')
args <- commandArgs(TRUE)
#exposure_name=args[1]
#outcome_name=args[2]
exposure_name="OverallAccelerationAverage"
outcome_name="CvdStatus"


#NOTE: TO MATCH THE GWAS 
covs=read.table("mr_covariates.euro.txt",header=TRUE,row.names=1,sep='\t')
covs$IID=NULL
covs[covs==-1000]=NA
covs$Sex=factor(covs$Sex)
covs$f.batch=factor(covs$f.batch)

snps=read.table(paste("mr_snp_subsets/",exposure_name,sep=""),header=TRUE,sep='\t',row.names=1)
snps[snps==-1]=NA

exposure_mr=read.table("mr_exposures.euro.txt",header=TRUE,row.names=1,sep='\t')
exposure_mr[exposure_mr==-1000]=NA
expo_v=subset(exposure_mr,select=c(exposure_name))

outcome_mr=read.table("mr_outcome.euro.txt",header=TRUE,row.names=1,sep='\t')
outcome_v=subset(outcome_mr,select=c(outcome_name))
outcome_v$CvdStatus=factor(outcome_v$CvdStatus)

x=subset(snps,select=c("rs59499656"))
y=expo_v


names(x) = c("x")
names(y) = c("y")
d = data.frame(x, y, covs)
#check whether the y is a factor (logistic regression) or a continuous value (linear regression)
if (is.factor(y) == TRUE)
{
  o = summary(glm(y ~ ., data = d, family = binomial))$coefficients
}else
{
  o = summary(lm(y ~ ., data = d))$coefficients
}
b = o["x", 1]
s = o["x", 2]
p = o["x", 4]

