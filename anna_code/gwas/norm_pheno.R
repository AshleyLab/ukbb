rm(list=ls())
library(ggplot2)
library(preprocessCore)

#read in matrix of phenotype traits; first 2 columns are FID & IID 
raw_pheno=read.table(" /scratch/PI/euan/projects/ukbb/gwas/accelerometry_plink/v1/accelerometery_aggregate_phenotypes.continuous.filtered.txt",header=TRUE,sep='\t')
cols=ncol(raw_pheno)
rows=nrow(raw_pheno)

#replace -1000 with NA for correct ranking behavior 
raw_pheno[raw_pheno==-1000]=NA

#remove obvious outliers for acceleration 
accel_outliers=which(raw_pheno[,3]>2000)
raw_pheno=raw_pheno[-c(accel_outliers),]
pheno=raw_pheno
for(i in seq(3,cols))
{
  print(i)
  pheno[,i]=(rank(pheno[,i],na.last=TRUE)-0.5)/rows
}

#perform quantile normalization 
pheno[,3:cols]=normalize.quantiles(as.matrix(pheno[,3:cols]))

#plot the raw vs normalized features 
library(ggplot2)
ggplot(data=pheno,aes(x=raw_pheno$OverallAccelerationAverage,y=pheno$OverallAccelerationAverage))+
  geom_point()+
  xlab("Overall Acceleration Average (mg)")+
  ylab("quantile.norm((rank - 0.5)/\nnumber_of_subjects)")+
  theme_bw(20)

#What if we just do quantile normalization on the values? 
pheno2=raw_pheno
pheno2[,3:cols]=normalize.quantiles(as.matrix(pheno2[,3:cols]))
ggplot(data=pheno,aes(x=raw_pheno$OverallAccelerationAverage,y=pheno2$OverallAccelerationAverage))+
  geom_point()+
  xlab("Overall Acceleration Average (mg)")+
  ylab("quantile.norm(Overall Acceleration Average)")+
  theme_bw(20)

#generate plots for all traits 
for(i in seq(3,cols))
{
  cur_trait=names(pheno)[i]
  png(filename=paste(cur_trait,'.png'))
  ggplot()
}