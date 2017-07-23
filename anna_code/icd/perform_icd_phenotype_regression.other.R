args = commandArgs(trailingOnly=TRUE)
col=as.numeric(args[1])
print(col) 
other=read.table("other.tsv",header=TRUE,sep='\t')
icd=read.table("icd.other.tsv",header=TRUE,sep='\t') 
cur_other=other[,col]
other_cols=ncol(other) 
icd_cols=ncol(icd) 

correction=152500
thresh=0.01 
significant_hits=data.frame(cur_other=character(),
		cur_icd=character(),
		pval=numeric(),
		odds=numeric())

for(c2 in seq(2,icd_cols))
{
cur_icd=factor(icd[,c2])
mod=glm(cur_icd ~ cur_other,family="binomial")
pval=(coef(summary(mod))[,4][2])*correction 
if (pval<thresh)
{
odds=exp(coef(summary(mod))[,1][2])
cur_other_name=names(other)[col] 
cur_icd_name=names(icd)[c2] 
cur_entry=data.frame(cur_other_name,cur_icd_name,pval,odds) 
significant_hits=rbind(significant_hits,cur_entry)
}
}

write.table(significant_hits,file=paste("significant_other_associations.tsv",args[1],sep=""),row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
