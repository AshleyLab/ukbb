rm(list=ls())
library(ggplot2)
library(preprocessCore)

#read in matrix of phenotype traits; first 2 columns are FID & IID 
raw_pheno=read.table("accelerometery_aggregate_phenotypes.continuous.filtered.txt",header=TRUE,sep='\t')
cols=ncol(raw_pheno)
rows=nrow(raw_pheno)

#replace -1000 with NA for correct ranking behavior 
raw_pheno[raw_pheno==-1000]=NA

#remove outliers more than 4 std dev away from the mean 
num_sd_permitted=4
outliers=c() 
for (i in seq(3,cols))
{
  mean_col=mean(na.omit(raw_pheno[,i]))
  sd_col=sd(na.omit(raw_pheno[,i]))
  upper_bound=mean_col+num_sd_permitted*sd_col 
  lower_bound=mean_col-num_sd_permitted*sd_col 
  outliers_high=which(raw_pheno[,i]>upper_bound)
  outliers_low=which(raw_pheno[,i]<lower_bound)
  outliers=c(outliers,outliers_high)
  outliers=c(outliers,outliers_low)
}
outliers=unique(sort(outliers))
raw_pheno=raw_pheno[-outliers,]

#commenting this out -- quantile normalization on the ranks seems to be giving an odd distribution; going to quantile noralize on values. 
#pheno=raw_pheno
# 
# for(i in seq(3,cols))
# {
#   pheno[,i]=(rank(pheno[,i],na.last=TRUE)-0.5)/rows
# }
# 
# #perform quantile normalization 
# pheno[,3:cols]=normalize.quantiles(as.matrix(pheno[,3:cols]))
# 
# #plot the raw vs normalized features 
# ggplot(data=pheno,aes(x=raw_pheno$OverallAccelerationAverage,y=pheno$OverallAccelerationAverage))+
#   geom_point()+
#   xlab("Overall Acceleration Average (mg)")+
#   ylab("quantile.norm((rank - 0.5)/\nnumber_of_subjects)")+
#   theme_bw(20)

#What if we just do quantile normalization on the values? 
pheno2=raw_pheno
pheno2[,3:cols]=normalize.quantiles(as.matrix(pheno2[,3:cols]))
pheno2=round(pheno2,2)
#write to output table 
pheno3=pheno2
pheno3[is.na(pheno3)]=-1000
write.table(pheno2,file="accelerometery_continuous_features_quantile_normalized_no_outliers.tsv",row.names = FALSE,col.names=TRUE,quote=FALSE,sep='\t')
browser() 
#sanity check -- scatter plots of original data w/ outliers removed vs. qnormed data 
for(i in seq(3,cols))
{
  cur_trait=names(pheno2)[i]
  png(filename=paste(cur_trait,'.png'))
  smoothScatter(raw_pheno[,i],pheno2[,i],xlab="Original",ylab="Quantile Normalized",main=cur_trait)
  dev.off() 
}
