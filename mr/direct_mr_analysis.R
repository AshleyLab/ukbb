# The flow from the RHR paper:
# 1. Rerun GWAS with p<1e-5, harsh LD threshold to increase power and independence
# 2. Used a single-sample MR

# fitness
# geno_paths = c('/Users/David/Desktop/ukbb/gwas_interpretation/mr/genotypes/fitness/oct_16_2017/') # NEW: Oct 2017
geno_paths = c('/Users/David/Desktop/ukbb/gwas_interpretation/mr/genotypes/fitness/oct_31_2017_exercise/',
               '/Users/David/Desktop/ukbb/gwas_interpretation/mr/genotypes/fitness/oct_31_2017_recovery/')# NEW: Oct 2017
pheno_path = '/Users/David/Desktop/ukbb/all_traits_for_gwas.RData'
# activity
#geno_path = '/Users/David/Desktop/ukbb/mr/genotypes/activity/'
#pheno_path = '/Users/David/Desktop/ukbb/fitness_analysis_final_fitness_scores.RData'
# other data
icd_path = '/Users/David/Desktop/ukbb/gwas_interpretation/mr/icd_matrix.txt'
snps_fitness_table = '/Users/David/Desktop/ukbb/top.snps.xlsx'
pcs_path = '/Users/David/Desktop/ukbb/covariates.augmented.txt'
covariates_path = '/Users/David/Desktop/ukbb/covariate_matrix.RData'
euro_ids_path = '/Users/David/Desktop/ukbb/pca_results_v2_chrom1_euro.eigenvec'

# # Read the genotypes from the raw files
# # ASSUMPTION: in each dir in geno_paths there are the same raw file names
# # with the same subjects, but different SNPs
# geno_path = geno_paths[1]
# geno_files = list.files(geno_path)
# geno_data = NULL
# for (geno_path in geno_paths){
#   curr_geno_files = list.files(geno_path)
#   curr_geno_files = curr_geno_files[grepl("\\.raw$",curr_geno_files)] 
#   curr_geno_data = NULL
#   for(f in curr_geno_files){
#     print(f)
#     snp_data = read.delim(paste(geno_path,f,sep=''),header=T,sep=" ")
#     rnames = snp_data[,1]
#     snp_data = snp_data[,-c(1:6)]
#     snp_names = colnames(snp_data)
#     inds = !grepl("_HET$",colnames(snp_data))
#     snp_data = snp_data[,inds]
#     snp_names = snp_names[inds]
#     if(length(snp_data)>0 && is.null(dim(snp_data))){snp_data = t(t(snp_data))}
#     rownames(snp_data) = rnames
#     colnames(snp_data) = snp_names
#     if(!is.null(curr_geno_data) && length(setdiff(rownames(curr_geno_data),rownames(snp_data)))>0){
#       print (paste("ERROR: row names do not match in: ", f))
#       next
#     }
#     if(is.null(curr_geno_data)){
#       curr_geno_data = snp_data
#       next
#     }
#     curr_geno_data = cbind(curr_geno_data,snp_data)
#   }
#   colnames(curr_geno_data) = gsub(colnames(curr_geno_data),pattern = "_.$",replace="")
#   print(dim(curr_geno_data))
#   if(!is.null(geno_data) && any(rownames(curr_geno_data)!=rownames(geno_data))){
#     print ("ERROR: row names do not match")
#     break
#   }
#   if(is.null(geno_data)){
#     geno_data = curr_geno_data
#   }
#   else{
#     geno_data = cbind(geno_data,curr_geno_data[,setdiff(colnames(curr_geno_data),colnames(geno_data))])
#   }
# }
# save(geno_data,file=paste(geno_path,"geno_data_from_raw_files.RData",sep=''))

# Save some time by loading directly:
geno_path = geno_paths[length(geno_paths)]
load(paste(geno_path,"geno_data_from_raw_files.RData",sep=''))
geno_data = as.matrix(geno_data)
mafs = apply(geno_data>0,2,function(x){x=x[!is.na(x)];mean(as.numeric(x>0))})
mafs = pmin(mafs,1-mafs)
hist(mafs,breaks=50)

gaus_norm<-function(x){
  x_r = (rank(x)-0.5)/length(x)
  x_n = qnorm(x_r)
  return(x_n)
}
IS_ACTIVITY = grepl(geno_path,pattern="activity")

# read the ids of the samples for the analysis: european
euro_pcs = read.delim(euro_ids_path)
euro_ids = as.character(euro_pcs$IID)
rm(euro_pcs);gc()

# read the exposure data
load(pheno_path)
if(IS_ACTIVITY){
  sort(apply(!is.na(accelerometry_scores),2,sum),decreasing=T)[1:10]
  expo = as.numeric(accelerometry_scores$NumberOfDaysModeratePhysicalActivity)
  names(expo) = rownames(accelerometry_scores)
  expo = expo[!is.na(expo)]
}
if(!IS_ACTIVITY){
  expo = fitness_scores_matrix
  # Use the original results
  # fuma_path = "/Users/David/Desktop/ukbb/gwas_interpretation/fuma/filtered_maf_0.01_1000genomes/pooled_p1_5e8_p2_1e4_ld_0.6_with_maf/"
  # exercise_snps = read.delim(paste(fuma_path,"leadSNPs.txt",sep=''),stringsAsFactors = F)$rsID
  # Use the new ones
  fuma_path1 = "/Users/David/Desktop/ukbb/gwas_interpretation/fuma/filtered_maf_0.01/pooled_exercise_hr_p1_5e6_ld_0.1_for_MR/"
  fuma_path2 = "/Users/David/Desktop/ukbb/gwas_interpretation/fuma/filtered_maf_0.01/pooled_recovery_p1_5e6_ld_0.1_for_MR/"
  snp_lists = list()
  snp_lists[["ExerciseHR"]] = read.delim(paste(fuma_path1,"leadSNPs.txt",sep=''),stringsAsFactors = F)$rsID
  snp_lists[["Recovery"]] = read.delim(paste(fuma_path2,"leadSNPs.txt",sep=''),stringsAsFactors = F)$rsID
  all_snps = unique(unlist(snp_lists))
  geno_data = geno_data[,intersect(all_snps,colnames(geno_data))]
  rm(accelerometry_scores)
}
# setdiff(unique(unlist(snp_lists)),colnames(geno_data))
# intersect(unique(unlist(snp_lists)),colnames(geno_data))
# setdiff(colnames(geno_data),unique(unlist(snp_lists)))

# read the clinical outcomes and partition the code into classes by:
# A00-B99  Certain infectious and parasitic diseases
# C00-D49  Neoplasms
# D50-D89  Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism
# E00-E89  Endocrine, nutritional and metabolic diseases
# F01-F99  Mental, Behavioral and Neurodevelopmental disorders
# G00-G99  Diseases of the nervous system
# H00-H59  Diseases of the eye and adnexa
# H60-H95  Diseases of the ear and mastoid process
# I00-I99  Diseases of the circulatory system
# J00-J99  Diseases of the respiratory system
# K00-K95  Diseases of the digestive system
# L00-L99  Diseases of the skin and subcutaneous tissue
# M00-M99  Diseases of the musculoskeletal system and connective tissue
# N00-N99  Diseases of the genitourinary system
# O00-O9A  Pregnancy, childbirth and the puerperium
# P00-P96  Certain conditions originating in the perinatal period
# Q00-Q99  Congenital malformations, deformations and chromosomal abnormalities
# R00-R99  Symptoms, signs and abnormal clinical and laboratory findings, not elsewhere classified
# S00-T88  Injury, poisoning and certain other consequences of external causes
# V00-Y99  External causes of morbidity
# Z00-Z99  Factors influencing health status and contact with health services
# Additional info from: http://www.healthnetworksolutions.net/index.php/understanding-the-icd-10-code-structure
# A & B:  Infectious and Parasitic Diseases
# C:  Neoplasms
# D:  Neoplasms, Blood, Blood-forming Organs
# E:  Endocrine, Nutritional, Metabolic
# F:  Mental and Behavioral Disorders
# G:  Nervous System
# H:  Eye and Adnexa, Ear and Mastoid Process
# I:  Circulatory System
# J:  Respiratory System
# K:  Digestive System
# L:  Skin and Subcutaneous Tissue
# M:  Musculoskeletal and Connective Tissue
# N:  Genitourinary System
# O:  Pregnancy, Childbirth and the Puerperium
# P:  Certain Conditions Originating in the Perinatal Period
# Q:  Congenital Malformations, Deformations and Chromosomal Abnormalities
# R:  Symptoms, Signs and Abnormal Clinical and Lab Findings
# S:  Injury, Poisoning, Certain Other Consequences of External Causes
# T:  Injury, Poisoning, Certain Other Consequences of External Causes
# U:  no codes listed, will be used for emergency code additions
# V, W, X, Y:  External Causes of Morbidity (homecare will only have to code how patient was hurt; other settings will also code where injury occurred, what activity patient was doing)
# Z:  Factors Influencing Health Status and Contact with Health Services (similar to current "V-codes")
icd_data = read.delim(icd_path)
rownames(icd_data) = as.character(icd_data[,1])
icd_data = icd_data[,-1]
icd_data_code_classes = list()
icd_data_code_classes[["heart_disease"]] = c(
  paste("I",21:22,sep=""), #myocardial infraction
  "I42","I50", #heart failure
  "I48", #atrial fibrilation
  "I471" #supraventricular tachycardia 
)
has_any_regex<-function(x,regexs){return(any(sapply(regexs,function(y,z)grepl(z,pattern=y),z=x)))}
icd_data_code_classes[["hypertension"]] = colnames(icd_data)[
  sapply(colnames(icd_data),has_any_regex,regexs = paste("^I",10:15,sep=""))]
icd_data_code_classes[["heart_disease"]] = colnames(icd_data)[
  sapply(colnames(icd_data),has_any_regex,regexs = c(paste("^I",21:22,sep=""),"I42","I50","I48","I471"))]
icd_data_code_classes[["diabetes"]] = colnames(icd_data)[
  sapply(colnames(icd_data),has_any_regex,regexs = paste("^E",11:13,sep=""))]
icd_data_code_classes[["cancer"]] =colnames(icd_data)[
  sapply(colnames(icd_data),nchar)==3 & grepl(colnames(icd_data),pattern="^C")
  ]
icd_data_code_classes[["disease_merge"]] = colnames(icd_data)[
  sapply(colnames(icd_data),nchar)==3 & grepl(colnames(icd_data),pattern="^(C|D|E|I|J|K|R)")
  ]
icd_data_code_columns = lapply(icd_data_code_classes,intersect,y=colnames(icd_data))
icd_data_code_columns = icd_data_code_columns[sapply(icd_data_code_columns,length)>0]
icd_outcome_matrices = lapply(icd_data_code_columns,function(x,y)y[,x],y=icd_data)
sapply(icd_outcome_matrices,dim)

healthy_subjs = rownames(icd_data)[rowSums(icd_data)==0]

# Read the covariates for the analysis
external_covs = read.delim(pcs_path)
table(external_covs$PC1!=(-1000))
# exclude samples with NAs in their PCs: these do not have genotypes
external_covs = external_covs[external_covs$PC1!=(-1000),]
pcs_matrix = external_covs[,grepl("^PC",colnames(external_covs))]
rownames(pcs_matrix) = external_covs$IID
batch_info = as.character(external_covs$f.batch)
batch_info[batch_info==-1000] = NA
pcs_matrix[pcs_matrix==-1000] = NA
names(batch_info) = external_covs$IID

# Merge the covariates and the pcs
load(covariates_path)
names_inters = intersect(rownames(pcs_matrix),rownames(covariate_matrix))
covariate_matrix = cbind(covariate_matrix[names_inters,],pcs_matrix[names_inters,],batch_info[names_inters])
colnames(covariate_matrix)[ncol(covariate_matrix)] = "batch"
pc_names = colnames(pcs_matrix)
simple_covs = covariate_matrix[,c("Sex","Age when attended assessment centre",pc_names,"batch")]
simple_covs[,"Sex"] = as.numeric(as.factor(simple_covs[,"Sex"]))-1
rm(external_covs);rm(pcs_matrix);gc()

# Define the subject sets and the data for the two-sample analysis
expo_subjs = rownames(expo)
if(is.null(expo_subjs)){expo_subjs = names(expo)}
outcome_subjs = setdiff(rownames(geno_data),expo_subjs)
one_sample_outcome_subjs = intersect(rownames(geno_data),expo_subjs)
subjects_with_data = intersect(rownames(simple_covs),euro_ids) # European that have covariates
subjects_with_data = intersect(subjects_with_data,rownames(geno_data)) # Have genotypes
expo_subjs = intersect(expo_subjs,subjects_with_data)
outcome_subjs = intersect(outcome_subjs,subjects_with_data)
outcome_subjs = intersect(outcome_subjs,rownames(icd_data))
one_sample_outcome_subjs = intersect(one_sample_outcome_subjs,subjects_with_data)
one_sample_outcome_subjs = intersect(one_sample_outcome_subjs,rownames(icd_data))
gc()
# Look at age and sex distribution
expo_subjs_age = simple_covs[expo_subjs,2]
nonexpo_subjs_age = simple_covs[setdiff(rownames(simple_covs),expo_subjs),2]
quantile(expo_subjs_age)
quantile(nonexpo_subjs_age)

# Define the tested outcomes
get_outcome_v<-function(x,y){
  x = as.matrix(x[intersect(rownames(x),y),])
  x[is.na(x)] = 0
  return(rowSums(x))
}
transform_outcome_to_binary<-function(x){
  nns = names(x)
  x = as.numeric(x>0)
  names(x) = nns
  return(x)
}
outcome_vs = lapply(icd_outcome_matrices,get_outcome_v,y=outcome_subjs)
par(mfrow=c(2,3));sapply(outcome_vs,hist);names(outcome_vs)
outcome_vs[["rhr"]] = gaus_norm(additional_scores$`Pulse rate`)
outcome_vs[["rhr"]] = outcome_vs[["rhr"]][intersect(names(outcome_vs[["rhr"]]),outcome_subjs)]
outcome_vs[1:4] = lapply(outcome_vs[1:4],transform_outcome_to_binary)
sapply(outcome_vs,length)
one_sample_outcome_vs = lapply(icd_outcome_matrices,get_outcome_v,y=one_sample_outcome_subjs)
par(mfrow=c(2,3));sapply(one_sample_outcome_vs,hist);names(one_sample_outcome_vs)
one_sample_outcome_vs[["rhr"]] = gaus_norm(additional_scores$`Pulse rate`)
one_sample_outcome_vs[["rhr"]] = one_sample_outcome_vs[["rhr"]][intersect(names(one_sample_outcome_vs[["rhr"]]),one_sample_outcome_subjs)]
one_sample_outcome_vs[1:4] = lapply(one_sample_outcome_vs[1:4],transform_outcome_to_binary)
sapply(one_sample_outcome_vs,length)
sapply(one_sample_outcome_vs[1:4],function(x)sum(x)/length(x))
sapply(outcome_vs[1:4],function(x)sum(x)/length(x))

check_maf_issue<-function(x){
  if((sum(x,na.rm=T)/2*sum(!is.na(x)))>0.5){
    print("fixing")
    return(2-x)
  }
  return(x)
}
corrected_maf_geno_data = apply(geno_data,2,check_maf_issue)
table(geno_data[,1],corrected_maf_geno_data[,1])

# # Oct 2017
# # Before MR, check if any of our SNPs is associated with selected traits
# # such as hemoglobin
# load(covariates_path)
# y1 = covariate_matrix[expo_subjs,"Haemoglobin concentration"]
# names(y1) = expo_subjs
# x1 = simple_covs[expo_subjs,]
# d = data.frame(y1,x1)
# d$batch = factor(d$batch)
# x1_lm = lm(y1~.,data=d)
# summary(x1_lm)
# hemo_resids = x1_lm$residuals
# plot(hemo_resids,y1[names(hemo_resids)])
# pvals = c()
# for(j in 1:ncol(geno_data)){
#   x2 = geno_data[names(hemo_resids),j]
#   curr_lm = lm(hemo_resids~x2)
#   curr_sum = summary(curr_lm)
#   p = curr_sum$coefficients[2,4]
#   pvals[j] = p
# }
# names(pvals) = colnames(geno_data)
# pvals[pvals<1e-4]
# mafs["rs16891982"]

# # Simple analyses of our scores vs the disease states
# icd_data_euro = icd_data[intersect(rownames(icd_data),subjects_with_data),]
# chisq_tests = list()
# for (nn in unique(unlist(icd_data_code_columns))){
#   ks_tests[[nn]] = list()
#   for(ss in colnames(geno_data)){
#     x1 = icd_data_euro[,nn]
#     x2 = geno_data[rownames(icd_data_euro),ss]
#     chisq_tests[[nn]][[ss]] = chisq.test(table(x1,x2))
#   }
# }
# chisq_pvals = lapply(chisq_tests,function(x)sapply(x,function(y)y$p.value))
# sapply(chisq_pvals,function(x)sum(x<0.01))

# Read the snp lists
if(!IS_ACTIVITY){
  # OLD CODE: keot to reproduce the old results in case we want to recheck
  # library(xlsx)
  # # dary'ls old summary table
  # fitness_snp_data = as.matrix(read.xlsx2(snps_fitness_table,sheetIndex = 1))
  # snp_lists = list()
  # snp_lists[['fitness_hr']] = unique(fitness_snp_data[fitness_snp_data[,"ancestry"]=="european" &
  #     (grepl(fitness_snp_data[,"trait"],pattern="Predicted_HR") | grepl(fitness_snp_data[,"trait"],pattern="slopes") | grepl(fitness_snp_data[,"trait"],pattern="Max_ach")) & 
  #     fitness_snp_data[,"score_adjustment"] == "simple"
  #     ,"snp"])
  # snp_lists[['exercise']] = unique(fitness_snp_data[fitness_snp_data[,"ancestry"]=="european" &
  #     (grepl(fitness_snp_data[,"trait"],pattern="Predicted_HR")) & 
  #     fitness_snp_data[,"score_adjustment"] == "simple"
  #     ,"snp"])
  # snp_lists[['recovery']] = fitness_snp_data[fitness_snp_data[,"ancestry"]=="european" &
  #     (grepl(fitness_snp_data[,"trait"],pattern="Rest")) & 
  #     fitness_snp_data[,"score_adjustment"] == "simple"
  #     ,"snp"]
  # # New runs from august 2017
  # snps_fitness_table = '/Users/David/Desktop/ukbb/gwas/daryl_clump/simple_norm_top.snps.xlsx'
  # fitness_snps = unique(as.character(read.xlsx2(snps_fitness_table,sheetIndex = 1)$snp))
  # snps_fitness_table = '/Users/David/Desktop/ukbb/gwas/daryl_clump/conservative_norm_top.snps.xlsx'
  # fitness_snps = unique(as.character(read.xlsx2(snps_fitness_table,sheetIndex = 1)$snp))
  # print(length(intersect(fitness_snps,colnames(geno_data))))
  
  # Fitness analysis: format all data to be on the same subjects
  expo_vs = list()
  expo_vs[["ExerciseHR"]] = expo[expo_subjs,2]
  expo_vs[["Recovery"]] = expo[expo_subjs,1]
  expo_vs = lapply(expo_vs,function(x)x[!is.na(x)])
  expo_vs = lapply(expo_vs,gaus_norm)
  
#   mr_analysis_res = list()
#   for(nn2 in names(expo_vs)){
#     for(nn3 in names(outcome_vs)){
#       currname = paste(nn2,nn3,"twosample",sep=",")
#       if(is.element(currname,set=names(mr_analysis_res))){next}
#       print(currname)
#       expo_v = expo_vs[[nn2]]
#       outcome_v = outcome_vs[[nn3]]
#       currsnps = intersect(snp_lists[[nn2]],colnames(geno_data))
#       mr_analysis_res[[currname]] = 
#           run_mr_analysis(expo_v,outcome_v,geno_data[,currsnps],covs=simple_covs,min_effect_size = 0)
#       print(mr_analysis_res[[currname]]$multivar_egger)
#       
#       currname = paste(nn2,nn3,"onesample",sep=",")
#       if(is.element(currname,set=names(mr_analysis_res))){next}
#       print(currname)
#       expo_v = expo_vs[[nn2]]
#       outcome_v = one_sample_outcome_vs[[nn3]]
#       inds = intersect(names(expo_v),names(outcome_v))
#       print(chisq.test(table(cut(expo_v[inds],5),outcome_v[inds]))$p.value)
#       currsnps = intersect(snp_lists[[nn2]],colnames(geno_data))
#       mr_analysis_res[[currname]] = 
#         run_mr_analysis(expo_v,outcome_v,geno_data[,currsnps],covs=simple_covs,min_effect_size = 0)
#       print(mr_analysis_res[[currname]]$multivar_egger)
#     }
#   }
}
# save(mr_analysis_res,file="/Users/David/Desktop/ukbb/gwas_interpretation/mr/results_nov_1_2017.RData")
# # Create a slim version of the results - remove the lm objects
# sapply(mr_analysis_res[[1]],class)
# for (i in 1:length(mr_analysis_res)){
#   mr_analysis_res[[i]]$unival_outcome_lm = NULL
#   mr_analysis_res[[i]]$univar_expo_lm = NULL
# }
# gc()
# save(mr_analysis_res,file="/Users/David/Desktop/ukbb/gwas_interpretation/mr/results_nov_1_2017_slim.RData")

load("/Users/David/Desktop/ukbb/gwas_interpretation/mr/results_nov_1_2017_slim.RData")

# Get a summary table
univar_ps = sapply(mr_analysis_res,function(x)attr(x[["univar_ivw"]],"Pval"))
univar_est = sapply(mr_analysis_res,function(x)attr(x[["univar_ivw"]],"Estimate"))
median_ps = sapply(mr_analysis_res,function(x)attr(x[["multivar_med"]],"Pvalue"))
median_est = sapply(mr_analysis_res,function(x)attr(x[["multivar_med"]],"Estimate"))
egger_ps = sapply(mr_analysis_res,function(x)attr(x[["multivar_egger"]],"Causal.pval"))
egger_est = sapply(mr_analysis_res,function(x)attr(x[["multivar_egger"]],"Estimate"))
mat = c()
GWS = paste(format(univar_est,digits = 2),paste('(',format(univar_ps,digits=3),')',sep=''),sep=' ')
Median = paste(format(median_est,digits = 2),paste('(',format(median_ps,digits=3),')',sep=''),sep=' ')
Egger = paste(format(egger_est,digits = 2),paste('(',format(egger_ps,digits=3),')',sep=''),sep=' ')
mat = cbind(GWS,Median,Egger)
rownames(mat) = names(univar_ps)

# Inspect a specific genetic weighted score: heart disease has inverted results
names(mr_analysis_res)
nn = names(mr_analysis_res)[1]
bx = mr_analysis_res[[nn]]$multivar_input@betaX
mm = as.matrix(geno_data[,names(bx)])
mm[is.na(mm)] = 0
gws = mm %*% bx
rhr = additional_scores$`Pulse rate`
names(gws) = rownames(geno_data)
out_name = strsplit(nn,split=',')[[1]][2]
inds = intersect(names(gws),names(outcome_vs[[out_name]]))
out_v = outcome_vs[[out_name]][inds];xx = gws[inds]

# For each disease look at sick vs healthy:
# Exercise HR
# RHR
# gws in exercise and non exercise subjects
icd_cols_corr_analysis = c()
get_disease_vs_scores_est<-function(x,y,sampsize=1000){
  inds = intersect(names(x),names(y))
  x = x[inds];y=y[inds]
  x1 = x[y==0]
  x2 = x[y==1]
  if(length(x2)<5){return(c(NA,NA,NA))}
  lmobj = summary(lm(x~y))$coefficients
  beta = lmobj[2,1]
  betap = lmobj[2,4]
  samp1 = sample(x1,min(length(x1),sampsize))
  samp2=sample(x2,min(length(x2),sampsize))
  pval = wilcox.test(samp1,samp2)$p.value
  return(c(beta,betap,pval))
}
library(speedglm)
get_disease_vs_scores_est_with_covs<-function(x,y,covs,useglm=F){
  inds = intersect(names(x),names(y))
  inds = intersect(inds,rownames(covs))
  x = x[inds];y=y[inds];covs=covs[inds,]
  d=data.frame(x=x,y=y,covs=covs)
  if(!useglm){
    lmobj = summary(lm(x~.,data=d))$coefficients
    betay = lmobj["y",1]
    betayp = lmobj["y",4]
    return(c(betay,betayp))
  }
  logiobj = speedglm(factor(y)~.,family=binomial(link='logit'),data=d)
  logiobj = summary(logiobj)
  logiobj = as.matrix(logiobj$coefficients)
  mode(logiobj) = 'numeric'
  betax = logiobj["x",1]
  betaxp = logiobj["x",4]
  return(c(betax,betaxp))
}
for(j in 1:ncol(icd_data)){
  disease_v = icd_data[,j]
  names(disease_v) = rownames(icd_data)
  expo_res = get_disease_vs_scores_est(expo_vs$ExerciseHR,disease_v)
  gws_res = get_disease_vs_scores_est(gws,disease_v)
  rhr_res = get_disease_vs_scores_est(rhr,disease_v)
  icd_cols_corr_analysis = rbind(icd_cols_corr_analysis,
    c(expo_res,gws_res,rhr_res))
  rownames(icd_cols_corr_analysis)[nrow(icd_cols_corr_analysis)] = colnames(icd_data)[j]
  print(colnames(icd_data)[j])
  print(rhr_res)
}
icd_cols_corr_analysis["R001",]
ord = order(icd_cols_corr_analysis[,9])
icd_cols_corr_analysis[ord,]
ord = order(icd_cols_corr_analysis[,6])
icd_cols_corr_analysis[ord,]

icd_cols_adjusted_corr_analysis = c()
batch_inds = grepl(colnames(simple_covs),pattern="batch")
for(j in 1:ncol(icd_data)){
  disease_v = icd_data[,j]
  names(disease_v) = rownames(icd_data)
  expo_res = get_disease_vs_scores_est_with_covs(expo_vs$ExerciseHR,disease_v,simple_covs[,!batch_inds])
  gws_res = get_disease_vs_scores_est_with_covs(gws,disease_v,simple_covs[,!batch_inds])
  rhr_res = get_disease_vs_scores_est_with_covs(rhr,disease_v,simple_covs[,!batch_inds])
  icd_cols_adjusted_corr_analysis = rbind(icd_cols_adjusted_corr_analysis,
                                 c(expo_res,gws_res,rhr_res))
  rownames(icd_cols_adjusted_corr_analysis)[nrow(icd_cols_adjusted_corr_analysis)] = colnames(icd_data)[j]
  print(colnames(icd_data)[j])
  print(rhr_res)
}

# Main features for exploratory analysis
brady = icd_data[,"R001"];names(brady)=rownames(icd_data)
expo_name = strsplit(nn,split=',')[[1]][1]
expo_v = expo_vs[[expo_name]]
outcome_v = transform_outcome_to_binary(rowSums(icd_outcome_matrices[[out_name]]))
bradyY = paste(outcome_v,brady[names(outcome_v)],sep=',')
names(bradyY) = names(outcome_v)
boxplot(rhr[names(outcome_v)]~bradyY)


out_v2 = one_sample_outcome_vs[[out_name]]
inds2 = intersect(names(out_v2),names(expo_v))
inds2 = names(expo_v)[!is.na(expo_v)]
G = gws[inds2];X=expo_v[inds2];Y=out_v2[inds2]
par(mfrow=c(1,3));boxplot(G~Y);boxplot(X~Y);boxplot(rhr[names(Y)]~Y)
table(Y,brady[inds2])
bradyY = paste(Y,brady[inds2],sep=',')
boxplot(X~bradyY)

load("/Users/david/Desktop/ukbb/subject2exercise_category.RData")
Z = subj2category[inds2]
cat1_subjects = inds2[grepl(pattern="Category 2",Z[inds2])]
boxplot(X[cat1_subjects]~Y[cat1_subjects])

M = icd_outcome_matrices[[out_name]][inds2,]
G2 = gws[setdiff(inds,inds2)];Y2=out_v[setdiff(inds,inds2)]

# Show disease vs. RHR in different age and sex strata
AGE = simple_covs[,2];names(AGE) = rownames(simple_covs)
AGE2 = gaus_norm(AGE)^2
U2 = simple_covs[names(YY),1:2]
U2[,2] = cut(U2[,2],5)
ZZ = paste("Sex:",U2[,1],", Age:",U2[,2],sep='')
table(ZZ)
AC1 = covariate_matrix[,"Number of days/week of moderate physical activity 10+ minutes"]
AC2 = covariate_matrix[,"Number of days/week of vigorous physical activity 10+ minutes"]
AC3 = covariate_matrix[,"Time spend outdoors in summer"]
AC4 = covariate_matrix[,"Number of days/week walked 10+ minutes"]
names(AC1) = rownames(covariate_matrix);names(AC2) = rownames(covariate_matrix);names(AC3) = rownames(covariate_matrix)
par(mfrow=c(2,5))
for(strat in sort(unique(ZZ))){
  inds = names(YY)[ZZ==strat]
  curr_rhr = rhr[inds]
  curr_y = YY[inds]
  curr_tt = table(curr_y)
  curr_names = paste(paste(names(curr_tt),",N=",curr_tt,sep=''))
  wilcox.test(AGE[inds][curr_y==1],AGE[inds][curr_y==0])$p.value
  boxplot(curr_rhr~curr_y,main=strat,ylim=c(30,140),
        names=curr_names,ylab="RHR")
}
for(strat in sort(unique(ZZ))){
  inds = names(YY)[ZZ==strat]
  curr_rhr = AC3[inds]
  inds = names(curr_rhr)[curr_rhr>=0]
  curr_rhr = curr_rhr[inds]
  curr_y = YY[inds]
  curr_tt = table(curr_y)
  curr_names = paste(paste(names(curr_tt),",N=",curr_tt,sep=''))
  boxplot(curr_rhr~curr_y,main=strat,names=curr_names,ylab="Activity")
}
inds = intersect(names(outcome_v),names(rhr))
d = data.frame(rhr=rhr[inds],Sex=simple_covs[inds,1],
               Age=simple_covs[inds,2],
               Age2 = AGE2[inds],HeartDisease=YY[inds],
               AC1 = AC1[inds],AC2=AC2[inds],AC3=AC3[inds])
d2 = data.frame(rhr=rhr[inds],covs=covariate_matrix[inds,feature_is_numeric],
                HeartDisease=YY[inds])
boxplot(rhr~Sex,data=d)
boxplot(rhr~HeartDisease,data=d)
summary(lm(rhr~.,data=d))
tempsum = summary(lm(rhr~.,data=d2))

# Redo MR, age<50
age_thr=55
age_selected_subjects = rownames(simple_covs)[simple_covs[,2]<age_thr]
age_gws = gws[age_selected_subjects]
lm1_names = intersect(names(Y2),age_selected_subjects)
lm1 = lm(Y2[lm1_names]~age_gws[lm1_names]+U2[lm1_names,1]+U2[lm1_names,2])
lm2_names = intersect(names(Y),age_selected_subjects)
lm2 = lm(X[lm2_names]~age_gws[lm2_names])

# Compare the disease outcomes
x1 = outcome_vs$hypertension
x2 = outcome_vs$disease_merge
table(x1,x2>0)

mr_plot(mr_analysis_res$`ExerciseHR,heart_disease,onesample`$multivar_input)
mr_plot(mr_analysis_res$`ExerciseHR,hypertension,onesample`$multivar_input)
mr_plot(mr_analysis_res$`ExerciseHR,rhr,onesample`$multivar_input)
mr_plot(mr_analysis_res$`ExerciseHR,rhr,twosample`$multivar_input)

# # Sanity check: compare recovery-rhr analysis to the two-sample indirect analysis
# # recovery analysis
# input_file = '/Users/David/Desktop/ukbb/mr/gwas_effects/recovery_snps_all_effects.txt'
# effects_data = read.delim(input_file)
# mr_data = get_twosample_mr_input(effects_data) # this is from the twosample_mr_analysis.R script
# mr_in_recovery_direct = mr_analysis_res$`Recovery,rhr`$multivar_input
# # exercise hr
# input_file = '/Users/David/Desktop/ukbb/mr/gwas_effects/exercise_snps_all_effects.txt'
# effects_data = read.delim(input_file)
# mr_data = get_twosample_mr_input(effects_data) # this is from the twosample_mr_analysis.R script
# mr_in_recovery_direct = mr_analysis_res$`ExerciseHR,rhr`$multivar_input
# direct_bx = mr_in_recovery_direct@betaX
# indirect_bx = mr_data[names(direct_bx),"bx"]
# plot(direct_bx,indirect_bx);abline(0,1)
# direct_bxse = mr_in_recovery_direct@betaXse
# indirect_bxse = mr_data[names(direct_bx),"bxse"]
# plot(direct_bxse,indirect_bxse);abline(0,1)
# direct_by = mr_in_recovery_direct@betaY
# indirect_by = mr_data[names(direct_bx),"by"]
# plot(direct_by,indirect_by);abline(0,1)
# direct_byse = mr_in_recovery_direct@betaYse
# indirect_byse = mr_data[names(direct_bx),"byse"]
# plot(direct_byse,indirect_byse);abline(0,1)
# direct_bx[mafs[names(direct_bx)]>0.5] = -direct_bx[mafs[names(direct_bx)]>0.5]
# direct_by[mafs[names(direct_by)]>0.5] = -direct_by[mafs[names(direct_by)]>0.5]
# plot(direct_bx,direct_by);abline(0,1)
# plot(indirect_bx,indirect_by);abline(0,1)
# # in case of disagreement:
# problematic_snps = names(which(direct_bx*indirect_bx < 0))
# direct_bx[problematic_snps]
# mafs[problematic_snps]
# colSums(geno_data[,problematic_snps]>0)
# snps = geno_data[,problematic_snps]

# if(IS_ACTIVITY){
#   # Activity analysis: format all data to be on the same subjects
#   subjs = names(expo)
#   sort(colSums(icd_data[subjs,]))
#   subjs = intersect(subjs,rownames(geno_data))
#   expo_v = expo[subjs]
#   outcome_v = heart_disease_outcome[subjs]
#   table(outcome_v)
#   # vs fitness
#   outcome_v = heart_disease_outcome[subjs]
#   fitness_v = fitness_scores_matrix[,2]
#   fitness_v = fitness_v[!is.na(fitness_v)]
#   subjs = intersect(subjs,names(fitness_v))
#   expo_v = expo[subjs]
#   outcome_v = heart_disease_outcome[subjs]
#   #####
#   # general disease
#   outcome_v = as.numeric(rowSums(icd_data[
#     subjs,sapply(colnames(icd_data),nchar)==3 & grepl(colnames(icd_data),pattern="^(A|B|C|D|E|F|G|I|J|K|M)")
#     ]))
#   #####
#   snps = geno_data[subjs,]
#   snps = as.matrix(snps)
#   snps = snps[,apply(snps>0,2,sum,na.rm=T)>100]
#   dim(snps)
#   gc()
#   # sample snps for a better running time
#   original_snps = snps
#   samp = sample(1:ncol(original_snps))[1:500]
#   snps = original_snps[,samp]
# }

get_lm_stats<-function(x,y){
  o = summary(lm(y~x))$coefficients
  b = o[-1,1]
  s = o[-1,2]
  return(c(b,s))
}
get_lm_stats_with_covs<-function(x,y,covs){
  if(is.null(covs)){return (get_lm_stats(x,y))}
  d = data.frame(x,y,covs)
  d$batch = factor(d$batch)
  o = summary(lm(y~.,data=d))$coefficients
  b = o["x",1]
  s = o["x",2]
  p = o["x",4]
  return(c(b,s,p))
}
gaus_norm<-function(x){
  x_r = (rank(x)-0.5)/length(x)
  x_n = qnorm(x_r)
  return(x_n)
}

library('MendelianRandomization')
run_mr_analysis<-function(expo_v,outcome_v,snps,covs=NULL,plot_mr_in=T,min_effect_size=0){
  covs = covs[rownames(snps),]
  
  expo_snps = snps[names(expo_v),]
  outcome_snps = snps[names(outcome_v),]
  expo_covs = covs[names(expo_v),]
  outcome_covs = covs[names(outcome_v),]
  
  lm_res = apply(expo_snps,2,get_lm_stats_with_covs,y=expo_v,covs=expo_covs)
  bx = lm_res[1,]
  bxse = lm_res[2,]
  lm_res = apply(outcome_snps,2,get_lm_stats_with_covs,y=outcome_v,covs=outcome_covs)
  by = lm_res[1,]
  byse = lm_res[2,]
  print("Two sample estimate completed")
  
  to_keep = abs(bx)>=0
  to_keep = abs(bx)>=min_effect_size
  if(sum(to_keep)==0){to_keep = abs(bx)>=median(abs(bx))}
  outcome_snps[is.na(outcome_snps)]=0
  expo_snps[is.na(expo_snps)]=0
  expo_snps = as.matrix(expo_snps)
  outcome_snps = as.matrix(outcome_snps)
  
  weighted_causal_v = expo_snps[,to_keep] %*% bx[to_keep]
  lm_expo_vs_g = lm(expo_v~weighted_causal_v)
  bx_obj = summary(lm_expo_vs_g)$coefficients
  bx_u = bx_obj[-1,1]
  bxse_u = bx_obj[-1,2]
  weighted_causal_v = outcome_snps[,to_keep] %*% bx[to_keep]
  lm_outcome_vs_g = lm(outcome_v~weighted_causal_v)
  by_obj = summary(lm_outcome_vs_g)$coefficients
  by_u = by_obj[-1,1]
  byse_u = by_obj[-1,2]
  
  res = list()
  if(sum(to_keep)>1){
    mr_in = mr_input(bx[to_keep],bxse[to_keep],by[to_keep],byse[to_keep])
    if(plot_mr_in){mr_plot(mr_in)}
    res[["multivar_input"]] = mr_in
    res[["multivar_egger"]] = mr_egger(mr_in,F,T)
    res[["multivar_med"]] = mr_median(mr_in)
    res[["multivar_all"]] = mr_allmethods(mr_in)
  }
  mr_in = mr_input(bx_u,bxse_u,by_u,byse_u)
  res[["univar_input"]] = mr_in
  res[["univar_ivw"]] = mr_ivw(mr_in)
  res[["univar_ml"]] = mr_maxlik(mr_in)
  res[["univar_expo_lm"]] = lm_expo_vs_g
  res[["unival_outcome_lm"]] = lm_outcome_vs_g
  return(res)
}

# Conditional independence tests
#install.packages("bnlearn")
library(bnlearn)
# d must have x,y, and z
run_discrete_ci_test<-function(d,cutsize=5){
  d = d[!apply(is.na(d),1,any),]
  for(nn in names(d)){
    if(is.numeric(d[[nn]])){
      print(nn)
      d[[nn]] = factor(cut(d[[nn]],breaks=cutsize))
    }
  }
  return(ci.test(d))
}

d = data.frame(x=weighted_causal_v,y=outcome_v,z=expo_v)
ci.test(d[!apply(is.na(d),1,any),])
run_discrete_ci_test(d)$p.value

d = data.frame(x=weighted_causal_v,y=expo_v,z=outcome_v)
ci.test(d[!apply(is.na(d),1,any),])$p.value
run_discrete_ci_test(d,2)$p.value
sum(!apply(is.na(d),1,any))

# # Get the data from the ped and map files (instead of raw)
# library(HardyWeinberg)
# transform_geno<-function(v,return_numeric=T){
#   v = as.matrix(v)
#   major = names(sort(table(v[,1]),decreasing=T))[1]
#   v[v==major]="A"
#   v[v!="A" & v!="0"]="B"
#   v[v=="0"] = "O"
#   ab = apply(v,1,paste,collapse='')
#   ab[ab=="BA"] = "AB"
#   #hw_test = HWLratio(table(ab)[c("AA","AB","BB")])
#   if(return_numeric){
#     ab[ab=="AA"] = "0"
#     ab[ab == "AB"] = "1"
#     ab[ab=="BB"] = "2"
#     ab[ab=="OO"] = NA
#     return(as.numeric(ab))
#   }
#   return(ab)
# }
# 
# # read the genotype data, get the current subjects at a time and merge the genotypes
# library(Matrix)
# geno_files = list.files(geno_path)
# geno_files = geno_files[!grepl(geno_files,pattern="RData$")]
# geno_names = unique(sapply(geno_files,function(x)paste(strsplit(x,split='\\.')[[1]][1:2],collapse='.')))
# geno_parsed_data = NULL
# for(g in geno_names){
#   print(g)
#   snp_data = read.delim(paste(geno_path,g,'.map',sep=''),header=F)
#   rdata = read.delim(paste(geno_path,g,'.ped',sep=''),sep=" ",header = F)
#   colnames(rdata) = c("famid", "pid", "fatid", "motid", "sex","affected",sapply(as.character(snp_data[,2]),rep,times=2))
#   if(ncol(rdata) != 6 + 2*nrow(snp_data)){
#     print ("ERROR: number of snps does not match the ped data")
#     break
#   }
#   gdata = c()
#   for (snp in as.character(snp_data[,2])){
#     inds = colnames(rdata)==snp
#     gv = transform_geno(rdata[,inds])
#     gdata = cbind(gdata,gv)
#     colnames(gdata)[ncol(gdata)] = snp
#   }
#   rownames(gdata) = rdata[,1]
#   gdata = Matrix(gdata)
#   if(!is.null(geno_parsed_data) && any(rownames(geno_parsed_data)!=rownames(gdata))){
#     print("ERROR:rownames do not match")
#     break
#   }
#   if(is.null(geno_parsed_data)){
#     geno_parsed_data = gdata
#   }
#   else{
#     geno_parsed_data = cbind(geno_parsed_data,gdata)
#   }
#   gc()
# }
# geno_data = geno_parsed_data
# save(geno_data,file=paste(geno_path,"fitness_geno_data_from_ped_files.RData",sep=''))
