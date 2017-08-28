# On Sherlock
geno_path = '/oak/stanford/groups/euan/projects/ukbb/gwas/snp_extract/fitness'
pheno_path = ''
icd_path = '/oak/stanford/groups/euan/projects/ukbb/code/anna_code/icd/icd_matrix.txt'

# Local
# fitness
geno_path = '/Users/David/Desktop/ukbb/mr/genotypes/fitness/'
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
  geno_file = '/Users/David/Desktop/ukbb/mr/genotypes/activity/geno_data_from_raw_files.RData'
}
if(!IS_ACTIVITY){
  geno_file = '/Users/David/Desktop/ukbb/mr/genotypes/fitness/geno_data_from_raw_files.RData'
  expo = fitness_scores_matrix
}

# read the genotype data
load(geno_file)
dim(geno_data)

# # read the genotypes from the raw files
# geno_files = list.files(geno_path)
# geno_files = geno_files[grepl("\\.raw$",geno_files)]
# geno_data = NULL
# for(f in geno_files){
#   print(f)
#   snp_data = read.delim(paste(geno_path,f,sep=''),header=T,sep=" ")
#   rownames(snp_data) = snp_data[,1]
#   snp_data = snp_data[,-c(1:6)]
#   snp_data = snp_data[,!grepl("_HET$",colnames(snp_data))]
#   if(!is.null(geno_data) && any(rownames(geno_data)!=rownames(snp_data))){
#     print ("ERROR: row names do not match")
#     break
#   }
#   if(is.null(geno_data)){
#     geno_data = snp_data
#     next
#   }
#   geno_data = cbind(geno_data,snp_data)
# }
# colnames(geno_data) = gsub(colnames(geno_data),pattern = "_.$",replace="")
# save(geno_data,file=paste(geno_path,"geno_data_from_raw_files.RData",sep=''))

# read the outcome
icd_data = read.delim(icd_path)
rownames(icd_data) = icd_data[,1]
icd_data = icd_data[,-1]
heart_disease_codes = c(
  paste("I",21:22,sep=""), #myocardial infraction
  "I42","I50", #heart failure
  "I48", #atrial fibrilation
  "I471" #supraventricular tachycardia 
)
hypertension_codes = paste("I",10:15,sep="")
diabetes_codes = paste("IE",10:14,sep="")
heart_disease_outcome = as.numeric(apply(icd_data[,intersect(colnames(icd_data),heart_disease_codes)],1,any))
hypertension_outcome = as.numeric(icd_data[,intersect(colnames(icd_data),hypertension_codes)])
names(heart_disease_outcome) = rownames(icd_data)
names(hypertension_outcome) = rownames(icd_data)

# read the covariates for the analysis
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

# merge the covariates and the pcs
load(covariates_path)
names_inters = intersect(rownames(pcs_matrix),rownames(covariate_matrix))
covariate_matrix = cbind(covariate_matrix[names_inters,],pcs_matrix[names_inters,],batch_info[names_inters])
colnames(covariate_matrix)[ncol(covariate_matrix)] = "batch"
pc_names = colnames(pcs_matrix)
simple_covs = covariate_matrix[,c("Sex","Age when attended assessment centre",pc_names,"batch")]
simple_covs[,"Sex"] = as.numeric(as.factor(simple_covs[,"Sex"]))-1
rm(external_covs);rm(pcs_matrix);rm(covariate_matrix);gc()

# get the subject set of the analysis and reshape the data
subjs = rownames(expo)
if(is.null(subjs)){
  subjs = names(expo)
}
subjs = intersect(subjs,euro_ids)
subjs = intersect(subjs,rownames(simple_covs))
subjs = intersect(subjs,rownames(icd_data))
subjs = intersect(subjs,rownames(geno_data))
geno_data = geno_data[subjs,]
icd_data = icd_data[subjs,]
simple_covs = simple_covs[subjs,]
gc()

# Read the snp lists
if(!IS_ACTIVITY){
  library(xlsx)
  fitness_snp_data = as.matrix(read.xlsx2(snps_fitness_table,sheetIndex = 1))
  fitness_snps = fitness_snp_data[fitness_snp_data[,"ancestry"]=="european" &
      (grepl(fitness_snp_data[,"trait"],pattern="Predicted_HR") | grepl(fitness_snp_data[,"trait"],pattern="Rest"))& 
      fitness_snp_data[,"score_adjustment"] == "conservative"
      ,"snp"]
  fitness_snps = fitness_snp_data[fitness_snp_data[,"ancestry"]=="european" &
      (grepl(fitness_snp_data[,"trait"],pattern="Predicted_HR") | grepl(fitness_snp_data[,"trait"],pattern="Rest")) & 
      fitness_snp_data[,"score_adjustment"] == "simple"
      ,"snp"]
  # Fitness analysis: format all data to be on the same subjects
  subjs = rownames(expo)[!is.na(expo[,2])]
  subjs = intersect(subjs,rownames(geno_data))
  expo_v = expo[subjs,2]
  outcome_v = heart_disease_outcome[subjs]
  # for sanity checks
  outcome_v = diabetes_outcome[subjs]
  # overal disease
  outcome_v = as.numeric(rowSums(icd_data[
      subjs,sapply(colnames(icd_data),nchar)==3 & grepl(colnames(icd_data),pattern="^(A|B|C|D|E|F|G|I|J|K|M)")
    ]))
  # vs rhr
  outcome_v = additional_scores$`Pulse rate`[subjs]
  table(outcome_v)
  snps = geno_data[subjs,intersect(fitness_snps,colnames(geno_data))]
  snps = as.matrix(snps)
  dim(snps)
  gc()
}

if(IS_ACTIVITY){
  # Activity analysis: format all data to be on the same subjects
  subjs = names(expo)
  sort(colSums(icd_data[subjs,]))
  subjs = intersect(subjs,rownames(geno_data))
  expo_v = expo[subjs]
  outcome_v = heart_disease_outcome[subjs]
  table(outcome_v)
  # vs fitness
  outcome_v = heart_disease_outcome[subjs]
  fitness_v = fitness_scores_matrix[,2]
  fitness_v = fitness_v[!is.na(fitness_v)]
  subjs = intersect(subjs,names(fitness_v))
  expo_v = expo[subjs]
  outcome_v = heart_disease_outcome[subjs]
  #####
  # general disease
  outcome_v = as.numeric(rowSums(icd_data[
    subjs,sapply(colnames(icd_data),nchar)==3 & grepl(colnames(icd_data),pattern="^(A|B|C|D|E|F|G|I|J|K|M)")
    ]))
  #####
  snps = geno_data[subjs,]
  snps = as.matrix(snps)
  snps = snps[,apply(snps>0,2,sum,na.rm=T)>100]
  dim(snps)
  gc()
  # sample snps for a better running time
  original_snps = snps
  samp = sample(1:ncol(original_snps))[1:500]
  snps = original_snps[,samp]
}

get_lm_stats<-function(x,y){
  o = summary(lm(y~x))$coefficients
  b = o[-1,1]
  s = o[-1,2]
  return(c(b,s))
}
get_lm_stats_with_covs<-function(x,y,covs){
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

expo_v = gaus_norm(expo_v)

# get betas and their sds, also merge snps based on betas
# # multivariate
# lm_expo_vs_g = lm(expo_v~snps)
# bx_obj = summary(lm_expo_vs_g)$coefficients
# bx = bx_obj[-1,1]
# bxse = bx_obj[-1,2]
# lm_outcome_vs_g = lm(outcome_v~snps)
# by_obj = summary(lm_outcome_vs_g)$coefficients
# by = by_obj[-1,1]
# byse = by_obj[-1,2]
# univariate
lm_res = apply(snps,2,get_lm_stats,y=expo_v)
bx = lm_res[1,]
bxse = lm_res[2,]
lm_res = apply(snps,2,get_lm_stats,y=outcome_v)
by = lm_res[1,]
byse = lm_res[2,]
# univariate with covs
covs = simple_covs[rownames(snps),]
lm_res = apply(snps,2,get_lm_stats_with_covs,y=expo_v,covs=covs)
bx = lm_res[1,]
bxse = lm_res[2,]
lm_res = apply(snps,2,get_lm_stats_with_covs,y=outcome_v,covs=covs)
by = lm_res[1,]
byse = lm_res[2,]

# # merge snps by beta and rerun
snps[is.na(snps)]=0
weighted_causal_v = snps %*% bx
hist(weighted_causal_v)
# lm_expo_vs_g = lm(expo_v~weighted_causal_v)
# bx_obj = summary(lm_expo_vs_g)$coefficients
# bx = bx_obj[-1,1]
# bxse = bx_obj[-1,2]
# lm_outcome_vs_g = lm(outcome_v~weighted_causal_v)
# by_obj = summary(lm_outcome_vs_g)$coefficients
# by = by_obj[-1,1]
# byse = by_obj[-1,2]
# 
# # outcome vs exposure
# o_v_e = glm(outcome_v~expo_v)
# summary(o_v_e)

# Perform MR
library('MendelianRandomization')
mr_in = mr_input(bx,bxse,by,byse)
mr_plot(mr_in)
res = list()
res[["maxlik"]] = mr_maxlik(mr_in)
res[["ivw"]] = mr_ivw(mr_in)
# multivariate analysis
res[["med"]] = mr_median(mr_in)
res[["egger"]] = mr_egger(mr_in)
res[["all"]] = mr_allmethods(mr_in)
attr(res$egger,"Pvalue.Est")
chisq.test(table(outcome_v,heart_disease_outcome[subjs]))$p.value

# Same tests on the top associated snps
num_top = 10
ord = order(abs(bx),decreasing = T)
mr_in = mr_input(bx[ord][1:num_top],bxse[ord][1:num_top],by[ord][1:num_top],byse[ord][1:num_top])
mr_plot(mr_in)
res = list()
res[["maxlik"]] = mr_maxlik(mr_in)
res[["ivw"]] = mr_ivw(mr_in)
# multivariate analysis
res[["med"]] = mr_median(mr_in)
res[["egger"]] = mr_egger(mr_in)
res[["all"]] = mr_allmethods(mr_in)
attr(res$egger,"Pvalue.Est")

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
