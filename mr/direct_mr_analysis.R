# On Sherlock
geno_path = '/oak/stanford/groups/euan/projects/ukbb/gwas/snp_extract/fitness'
pheno_path = ''
icd_path = '/oak/stanford/groups/euan/projects/ukbb/code/anna_code/icd/icd_matrix.txt'

# Local
# fitness
geno_path = '/Users/David/Desktop/ukbb/mr/genotypes/fitness/' # OLD
geno_path = '/Users/David/Desktop/ukbb/mr/genotypes/fitness/oct_16_2017/' # NEW: Oct 2017
pheno_path = '/Users/David/Desktop/ukbb/all_traits_for_gwas.RData'
# activity
#geno_path = '/Users/David/Desktop/ukbb/mr/genotypes/activity/'
#pheno_path = '/Users/David/Desktop/ukbb/fitness_analysis_final_fitness_scores.RData'
# other data
icd_path = '/Users/David/Desktop/ukbb/mr/icd_matrix.txt'
snps_fitness_table = '/Users/David/Desktop/ukbb/top.snps.xlsx'
pcs_path = '/Users/David/Desktop/ukbb/covariates.augmented.txt'
covariates_path = '/Users/David/Desktop/ukbb/covariate_matrix.RData'
euro_ids_path = '/Users/David/Desktop/ukbb/pca_results_v2_chrom1_euro.eigenvec'

# read the genotypes from the raw files
geno_files = list.files(geno_path)
geno_files = geno_files[grepl("\\.raw$",geno_files)]
geno_data = NULL
for(f in geno_files){
  print(f)
  snp_data = read.delim(paste(geno_path,f,sep=''),header=T,sep=" ")
  rnames = snp_data[,1]
  snp_data = snp_data[,-c(1:6)]
  snp_names = colnames(snp_data)
  inds = !grepl("_HET$",colnames(snp_data))
  snp_data = snp_data[,inds]
  snp_names = snp_names[inds]
  if(length(snp_data)>0 && is.null(dim(snp_data))){
    snp_data = t(t(snp_data))
  }
  rownames(snp_data) = rnames
  colnames(snp_data) = snp_names
  if(!is.null(geno_data) && any(rownames(geno_data)!=rownames(snp_data))){
    print ("ERROR: row names do not match")
    break
  }
  if(is.null(geno_data)){
    geno_data = snp_data
    next
  }
  geno_data = cbind(geno_data,snp_data)
}
colnames(geno_data) = gsub(colnames(geno_data),pattern = "_.$",replace="")
dim(geno_data)
save(geno_data,file=paste(geno_path,"geno_data_from_raw_files.RData",sep=''))

# Save some time by loading directly:
load(paste(geno_path,"geno_data_from_raw_files.RData",sep=''))
geno_data = as.matrix(geno_data)
mafs = apply(geno_data>0,2,function(x){x=x[!is.na(x)];mean(as.numeric(x>0))})

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
}

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
rownames(icd_data) = icd_data[,1]
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
rm(external_covs);rm(pcs_matrix);rm(covariate_matrix);gc()

# Define the subject sets and the data for the two-sample analysis
expo_subjs = rownames(expo)
if(is.null(expo_subjs)){expo_subjs = names(expo)}
outcome_subjs = setdiff(rownames(geno_data),expo_subjs)

subjects_with_data = intersect(rownames(simple_covs),euro_ids) # European that have covariates
subjects_with_data = intersect(subjects_with_data,rownames(geno_data)) # Have genotypes
expo_subjs = intersect(expo_subjs,subjects_with_data)
outcome_subjs = intersect(outcome_subjs,subjects_with_data)
outcome_subjs = intersect(outcome_subjs,rownames(icd_data))
gc()

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
outcome_vs[["rhr"]] = gaus_norm(additional_scores$`Pulse rate`)
outcome_vs[["rhr"]] = outcome_vs[["rhr"]][intersect(names(outcome_vs[["rhr"]]),outcome_subjs)]
outcome_vs[1:4] = lapply(outcome_vs[1:4],transform_outcome_to_binary)
sapply(outcome_vs,length)

# Oct 2017
# Before MR, check if any of our SNPs is associated with selected traits
# such as hemoglobin
load(covariates_path)
y1 = covariate_matrix[expo_subjs,"Haemoglobin concentration"]
names(y1) = expo_subjs
x1 = simple_covs[expo_subjs,]
d = data.frame(y1,x1)
d$batch = factor(d$batch)
x1_lm = lm(y1~.,data=d)
summary(x1_lm)
hemo_resids = x1_lm$residuals
plot(hemo_resids,y1[names(hemo_resids)])
pvals = c()
for(j in 1:ncol(geno_data)){
  x2 = geno_data[names(hemo_resids),j]
  curr_lm = lm(hemo_resids~x2)
  curr_sum = summary(curr_lm)
  p = curr_sum$coefficients[2,4]
  pvals[j] = p
}
names(pvals) = colnames(geno_data)
pvals[pvals<1e-4]

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
  
  fuma_path = "/Users/David/Desktop/ukbb/gwas_interpretation/fuma/pooled_p1_5e8_p2_1e4_ld_0.6_with_maf/"
  exercise_snps = read.delim(paste(fuma_path,"leadSNPs.txt",sep=''),stringsAsFactors = F)$rsID
  
  snp_lists = list()
  snp_lists[["ExerciseHR"]] = intersect(colnames(geno_data),exercise_snps)
  snp_lists[["Recovery"]] = setdiff(colnames(geno_data),exercise_snps)

  # Fitness analysis: format all data to be on the same subjects
  expo_vs = list()
  expo_vs[["ExerciseHR"]] = expo[expo_subjs,2]
  expo_vs[["Recovery"]] = expo[expo_subjs,1]
  expo_vs = lapply(expo_vs,function(x)x[!is.na(x)])
  expo_vs = lapply(expo_vs,gaus_norm)
  
  mr_analysis_res = list()
  for(nn2 in names(expo_vs)){
    for(nn3 in names(outcome_vs)){
      currname = paste(nn2,nn3,sep=",")
      if(is.element(currname,set=names(mr_analysis_res))){next}
      print(currname)
      expo_v = expo_vs[[nn2]]
      outcome_v = outcome_vs[[nn3]]
      currsnps = snp_lists[[nn2]]
      mr_analysis_res[[currname]] = 
          run_mr_analysis(expo_v,outcome_v,geno_data[,currsnps],covs=simple_covs,min_effect_size = 0)
      print(mr_analysis_res[[currname]]$multivar_egger)
    }
  }
}
p.adjust(sapply(mr_analysis_res,function(x)attr(x[["univar_ivw"]],"Pval")),method='fdr')
p.adjust(sapply(mr_analysis_res,function(x)attr(x[["multivar_egger"]],"Pvalue.Est")),method='fdr')

# Compare the disease outcomes
x1 = outcome_vs$hypertension
x2 = outcome_vs$disease_merge
table(x1,x2>0)

mr_plot(mr_analysis_res$`ExerciseHR,heart_disease`$multivar_input)
mr_plot(mr_analysis_res$`ExerciseHR,hypertension`$multivar_input)
mr_plot(mr_analysis_res$`ExerciseHR,rhr`$multivar_input)
save(mr_analysis_res,file="/Users/David/Desktop/ukbb/mr/results_oct_2017.RData")

# Sanity check: compare recovery-rhr analysis to the two-sample indirect analysis
# recovery analysis
input_file = '/Users/David/Desktop/ukbb/mr/gwas_effects/recovery_snps_all_effects.txt'
effects_data = read.delim(input_file)
mr_data = get_twosample_mr_input(effects_data) # this is from the twosample_mr_analysis.R script
mr_in_recovery_direct = mr_analysis_res$`Recovery,rhr`$multivar_input
# exercise hr
input_file = '/Users/David/Desktop/ukbb/mr/gwas_effects/exercise_snps_all_effects.txt'
effects_data = read.delim(input_file)
mr_data = get_twosample_mr_input(effects_data) # this is from the twosample_mr_analysis.R script
mr_in_recovery_direct = mr_analysis_res$`ExerciseHR,rhr`$multivar_input
direct_bx = mr_in_recovery_direct@betaX
indirect_bx = mr_data[names(direct_bx),"bx"]
plot(direct_bx,indirect_bx);abline(0,1)
direct_bxse = mr_in_recovery_direct@betaXse
indirect_bxse = mr_data[names(direct_bx),"bxse"]
plot(direct_bxse,indirect_bxse);abline(0,1)
direct_by = mr_in_recovery_direct@betaY
indirect_by = mr_data[names(direct_bx),"by"]
plot(direct_by,indirect_by);abline(0,1)
direct_byse = mr_in_recovery_direct@betaYse
indirect_byse = mr_data[names(direct_bx),"byse"]
plot(direct_byse,indirect_byse);abline(0,1)
direct_bx[mafs[names(direct_bx)]>0.5] = -direct_bx[mafs[names(direct_bx)]>0.5]
direct_by[mafs[names(direct_by)]>0.5] = -direct_by[mafs[names(direct_by)]>0.5]
plot(direct_bx,direct_by);abline(0,1)
plot(indirect_bx,indirect_by);abline(0,1)
# in case of disagreement:
problematic_snps = names(which(direct_bx*indirect_bx < 0))
direct_bx[problematic_snps]
mafs[problematic_snps]
colSums(geno_data[,problematic_snps]>0)
snps = geno_data[,problematic_snps]

names(mr_analysis_res)
name1 = "ExerciseHR,heart_disease"
mr_in = mr_analysis_res[[name1]][["multivar_input"]]
mr_plot(mr_in)
mr_analysis_res[[name1]][["multivar_egger"]]
mr_analysis_res[[name1]][["multivar_all"]]
expo_v = expo_vs[["ExerciseHR"]]
outcome_v = outcome_vs[["heart_disease"]]
subjs = intersect(names(expo_v),rownames(geno_data))
subjs = intersect(names(outcome_v),subjs)
expo_v = expo_v[subjs]
outcome_v = outcome_v[subjs]

d = data.frame(x=weighted_causal_v,y=outcome_v,z=expo_v)
run_discrete_ci_test(d,5)$p.value

d = data.frame(x=weighted_causal_v,y=expo_v,z=outcome_v)
run_discrete_ci_test(d,5)$p.value

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
  return(c(b,s))
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
  plot(bx,by)

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
    res[["multivar_egger"]] = mr_egger(mr_in,T,T)
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
run_discrete_ci_test<-function(d,cutsize=10){
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
