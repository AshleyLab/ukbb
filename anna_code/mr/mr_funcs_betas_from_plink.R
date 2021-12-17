rm(list=ls())
exposure_name="X6am_12pm"
outcome_name="CvdStatus"

#read in the PLINK results 
data=read.table("gwas.summary",header=TRUE,sep='\t',stringsAsFactors = FALSE)

#get the betas and se values from the PLINK regression 
exposure_subset=data[data$Trait==exposure_name,]
snps=exposure_subset$SNP 

outcome_subset= data[(data$Trait==outcome_name) & (data$SNP %in% snps) & (data$PVAL < .10),]
outcome_subset=outcome_subset[order(outcome_subset$SNP),]

snps=outcome_subset$SNP

exposure_subset=exposure_subset[(exposure_subset$SNP %in% snps),]
exposure_subset=exposure_subset[order(exposure_subset$SNP),]

exposure_beta=exposure_subset$BETA
exposure_se=exposure_subset$SE 
outcome_beta=outcome_subset$BETA 
outcome_se = outcome_subset$SE

#generate MendelianRandomization inputs 
mr_in = mr_input(exposure_beta, 
                 exposure_se, 
                 outcome_beta, 
                 outcome_se,
                 snps=snps,
                 exposure=exposure_name,
                 outcome=outcome_name)

#perform MR analysis 
MRAllObject_all <- mr_allmethods(mr_in, method = "all")
print(mr_plot(mr_in,interactive=FALSE,labels=TRUE))
print(mr_plot(MRAllObject_all,labels=TRUE))
MRAllObject_all

#if < 3 SNPs, only mr_ivw, mr_maxlik 