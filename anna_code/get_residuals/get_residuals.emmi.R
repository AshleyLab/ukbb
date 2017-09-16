rm(list=ls())
covar=data.frame(read.table('covariates.emmi.filtered.txt',header=TRUE,sep='\t'))
phenotypes_continuous=data.frame(read.table('ukb_PA_GS_SLEEP.pheno.recoded',header=TRUE,sep='\t'))


#replace -1000 with NA 
covar[covar==-1000]=NA
phenotypes_continuous[phenotypes_continuous==-1000]=NA


#sort all matrices by FID 
covar=covar[order(covar$FID),]
#remove subjects with negative id's 
fid_to_keep=covar$FID
iid_to_keep=covar$IID
covar$FID=NULL
covar$IID=NULL 
#identify factors & ordinals 
covar$Sex=factor(covar$Sex)
covar$f.batch=factor(covar$f.batch)

phenotypes_continuous=phenotypes_continuous[order(phenotypes_continuous$FID),]
pheno_names_continuous=names(phenotypes_continuous)
phenotypes_continuous$FID=NULL
phenotypes_continuous$IID=NULL

phenotypes_continuous=as.matrix(phenotypes_continuous)

residuals_continuous=matrix(nrow=nrow(phenotypes_continuous),ncol=ncol(phenotypes_continuous))
for(i in seq(1,ncol(phenotypes_continuous)))
{
curdata=cbind(phenotypes_continuous[,i],covar)
names(curdata)=c("Y",names(covar))
curdata$Sex=factor(curdata$Sex)
curdata$f.batch=factor(curdata$f.batch)
residuals_continuous[,i]=as.vector(residuals(lm(Y ~ Sex
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
             +f.batch,data=curdata,na.action=na.exclude),na.action=na.exclude))
print(i)
}

residuals_continuous=cbind(fid_to_keep,iid_to_keep,as.data.frame(residuals_continuous))
names(residuals_continuous)=pheno_names_continuous
write.table(x=round(residuals_continuous,2),file="ukb_PA_GS_SLEEP.pheno.recoded.residuals",quote=FALSE,row.names = FALSE,col.names = TRUE,sep = '\t')
