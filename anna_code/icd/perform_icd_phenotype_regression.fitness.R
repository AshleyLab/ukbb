args = commandArgs(trailingOnly=TRUE)
col=as.numeric(args[1])
print(col) 
fitness=read.table("fitness.tsv",header=TRUE,sep='\t')
icd=read.table("icd.fitness.tsv",header=TRUE,sep='\t') 
cur_fitness=fitness[,col]
fitness_cols=ncol(fitness) 
icd_cols=ncol(icd) 

correction=152500
thresh=0.01 
significant_hits=data.frame(cur_fitness=character(),
		cur_icd=character(),
		pval=numeric(),
		odds=numeric())

for(c2 in seq(2,icd_cols))
{
cur_icd=factor(icd[,c2])
mod=glm(cur_icd ~ cur_fitness,family="binomial")
pval=(coef(summary(mod))[,4][2])*correction 
if (pval<thresh)
{
odds=exp(coef(summary(mod))[,1][2])
cur_fitness_name=names(fitness)[col] 
cur_icd_name=names(icd)[c2] 
cur_entry=data.frame(cur_fitness_name,cur_icd_name,pval,odds) 
significant_hits=rbind(significant_hits,cur_entry)
}
}

write.table(significant_hits,file=paste("significant_fitness_associations.tsv",args[1],sep=""),row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
