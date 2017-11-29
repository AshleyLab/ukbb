# Functions
gaus_norm<-function(x){
  x_r = (rank(x)-0.5)/length(x)
  x_n = qnorm(x_r)
  return(x_n)
}
pairwise_cor<-function(x,y,...){
  inds = !is.na(x)&!is.na(y)
  return(cor(x[inds],y[inds],...))
}
get_lm_stats_with_covs<-function(x,y,covs){
  if(is.null(covs)){return (get_lm_stats(x,y))}
  d = data.frame(x,y,covs)
  if(is.element("batch",colnames(d))){d$batch = factor(d$batch)}
  o = summary(lm(y~.,data=d))$coefficients
  b = o["x",1]
  s = o["x",2]
  p = o["x",4]
  return(c(b,s,p))
}
analyze_source_vs_target_node<-function(source,target,addToS=NULL,data,ci_test=run_ci_test,depth=1,p0=1e-4,pthr=1e-4){
  pairwise_p = ci_test(source,target,NULL,data)
  if(pairwise_p>p0){return(list(res=F,sepset=c()))}
  if(depth<1){return(list(res=T,sepset=c()))}
  n = ncol(data)
  if(is.null(n)){n = ncol(data$DATA)}
  colinds = 1:n;colinds=colinds[-c(source,target,addToS)]
  n = length(colinds)
  for(setsize in 1:depth){
    S = 1:setsize;islast=F
    while(!islast){
      tmp = getNextSet(n,setsize,S)
      if(is.element(source,set=S)||is.element(target,set=S)){
        S = tmp$nextSet
        islast = tmp$wasLast
        next
      }
      p = ci_test(source,target,union(addToS,colinds[S]),data)
      if(p>pthr){return(list(res=F,sepset=union(addToS,colinds[S])))}
      S = tmp$nextSet
      islast = tmp$wasLast
    }
  }
  return(list(res=T,sepset=c()))
}
library(bnlearn)
# install.packages('pcalg')
library(pcalg);library(bnlearn)
disc_data_using_cut<-function(x,cuts=5,min_bin_size=100){
  y = NULL
  if(!is.numeric(x)){y = factor(x)}
  if(length(unique(x))<=cuts){y = factor(x)}
  if(is.null(y)){y=factor(cut(x,breaks=cuts, ordered_result=T))}
  table_y = table(y)
  while(any(table_y<min_bin_size) && !all(y==y[1],na.rm=T)){
    curr_levels = levels(y)
    j = which(table_y==min(table_y))[1]
    ll = names(j)
    j2 = j-1
    if(j==1){j2=2}
    ll2 = names(table_y)[j2]
    new_levels=curr_levels
    new_levels[j]=paste(ll2,ll,sep=",")
    new_levels[j2]=paste(ll2,ll,sep=",")
    levels(y) = new_levels
    table_y = table(y)
  }
  return(y)
}
run_discrete_ci_test<-function(x,y,z,data,test="mi-adf",...){
  if(is.numeric(x)){x = names(data)[x]}
  if(is.numeric(y)){y = names(data)[y]}
  if(is.numeric(z) && length(z)>0){
    zz = c()
    for(ii in z){zz = c(zz,names(data)[ii])}
    z = zz
  }
  if(length(z)==0){
    inds = !apply(is.na(data[,c(x,y)]),1,any)
    return(ci.test(x,y,data=data[inds,c(x,y)],test=test,...)$p.value)
  }
  inds = !apply(is.na(data[,c(x,y,z)]),1,any)
  return(ci.test(x,y,z,data[inds,c(x,y,z)],test=test,...)$p.value)
}
run_ci_test_one_is_numeric<-function(x,y,z,data){
  if(is.numeric(x)){x = names(data)[x]}
  if(is.numeric(y)){y = names(data)[y]}
  if(is.numeric(z) && length(z)>0){
    zz = c()
    for(ii in z){zz = c(zz,names(data)[ii])}
    z = zz
  }
  yv = data[,y];xv = data[,x];zzv = NULL
  if(length(z)>0){zzv = data[,z]}
  if(is.numeric(yv)&&is.numeric(xv)){
    if(length(z)==0){
      return(ci.test(x,y,data=data,test="cor")$p.value)
    }
    else{
      d1 = data.frame(x=xv,zzv);d2=data.frame(y=yv,zzv)
      lm1  = lm(x~.,data=d1)$residuals
      lm2  = lm(y~.,data=d2)$residuals
      return(cor.test(lm1,lm2)$p.value)
    }
  }
  if(is.numeric(yv)&&length(z)>0){
    summ = summary(lm(y~.,data=data.frame(y=yv,x=xv,zzv)))
    return(summ$coefficients[2,4])
  }
  if(is.numeric(xv)&&length(z)>0){
    summ = summary(lm(x~.,data=data.frame(x=xv,y=yv,zzv)))
    return(summ$coefficients[2,4])
  }
  if(is.numeric(yv)&&length(z)==0){
    summ = summary(lm(y~.,data=data.frame(y=yv,x=xv)))
    return(summ$coefficients[2,4])
  }
  if(is.numeric(xv)&&length(z)==0){
    summ = summary(lm(x~.,data=data.frame(x=xv,y=yv)))
    return(summ$coefficients[2,4])
  }
}
# data is a list with DATA and discDATA
run_ci_test<-function(x,y,z,data,test="mi-adf",...){
  yv = data$DATA[,y];xv = data$DATA[,x]
  if(is.numeric(yv)||is.numeric(xv)){
    return(run_ci_test_one_is_numeric(x,y,z,data$DATA))
  }
  return(run_discrete_ci_test(x,y,z,data$discDATA,test=test))
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

##################################################
##################################################
##################################################
# Load data
geno_paths = c('/Users/David/Desktop/ukbb/gwas_interpretation/mr/genotypes/fitness/oct_31_2017_exercise/',
               '/Users/David/Desktop/ukbb/gwas_interpretation/mr/genotypes/fitness/oct_31_2017_recovery/')# NEW: Oct 2017
pheno_path = '/Users/David/Desktop/ukbb/all_traits_for_gwas.RData'
icd_path = '/Users/David/Desktop/ukbb/gwas_interpretation/mr/icd_matrix.txt'
snps_fitness_table = '/Users/David/Desktop/ukbb/top.snps.xlsx'
pcs_path = '/Users/David/Desktop/ukbb/covariates.augmented.txt'
covariates_path = '/Users/David/Desktop/ukbb/covariate_matrix.RData'
euro_ids_path = '/Users/David/Desktop/ukbb/pca_results_v2_chrom1_euro.eigenvec'

# Load genotypes
geno_path = geno_paths[length(geno_paths)]
load(paste(geno_path,"geno_data_from_raw_files.RData",sep=''))
geno_data = as.matrix(geno_data)
# # transform data into mafs
# is_maf_represented = apply(geno_data,2,table)
# is_maf_represented = apply(is_maf_represented,2,function(x)x/sum(x))
# zero_allele_prop = is_maf_represented[1,] + 0.5*is_maf_represented[2,]
# min(zero_allele_prop)
# Impute missing values
library(impute)
im = impute.knn(geno_data,k = 10)

# read the ids of the samples for the analysis: european
euro_pcs = read.delim(euro_ids_path)
euro_ids = as.character(euro_pcs$IID)
rm(euro_pcs);gc()

# read the exposure data
load(pheno_path)
expo = fitness_scores_matrix
fuma_path1 = "/Users/David/Desktop/ukbb/gwas_interpretation/fuma/filtered_maf_0.01/pooled_exercise_hr_p1_5e6_ld_0.1_for_MR/"
fuma_path2 = "/Users/David/Desktop/ukbb/gwas_interpretation/fuma/filtered_maf_0.01/pooled_recovery_p1_5e6_ld_0.1_for_MR/"
snp_lists = list()
snp_lists[["ExerciseHR"]] = read.delim(paste(fuma_path1,"leadSNPs.txt",sep=''),stringsAsFactors = F)$rsID
snp_lists[["Recovery"]] = read.delim(paste(fuma_path2,"leadSNPs.txt",sep=''),stringsAsFactors = F)$rsID
all_snps = unique(unlist(snp_lists))
geno_data = geno_data[,intersect(all_snps,colnames(geno_data))]
rm(accelerometry_scores)
# Get gwascat data
gwascat_data = list()
gwascat_data[["ExerciseHR"]] = read.delim(paste(fuma_path1,"gwascatalog.txt",sep=''),stringsAsFactors = F)
gwascat_data[["Recovery"]] = read.delim(paste(fuma_path2,"gwascatalog.txt",sep=''),stringsAsFactors = F)
region_data = list()
region_data[["ExerciseHR"]] = read.delim(paste(fuma_path1,"leadSNPs.txt",sep=''),stringsAsFactors = F)
region_data[["Recovery"]] = read.delim(paste(fuma_path2,"leadSNPs.txt",sep=''),stringsAsFactors = F)
sort(table(gwascat_data[["ExerciseHR"]]$Trait))

# Get the phenotypes - all data
load("/Users/david/Desktop/ukbb/biobank_collated_pheno_data.RData")
medication_cols = colnames(pheno_data)[grepl(pattern="medica",colnames(pheno_data),ignore.case = T)]
medication_data = pheno_data[,medication_cols]
# Map the medications to beta-blockers and calcium channel blockers
bblockers = as.character(read.delim('/Users/David/Desktop/ukbb/list_of_beta_blockers.txt')[,1])
ccblockers = as.character(read.delim('/Users/David/Desktop/ukbb/list_of_calcium_channel_blockers.txt')[,1])
medcode2name = read.delim('/Users/David/Desktop/ukbb/treatment_coding4.tsv',row.names=1,stringsAsFactors = F)
ns = rownames(medcode2name);medcode2name=medcode2name[,1];names(medcode2name)=ns
name2medcode = ns;names(name2medcode)=unname(medcode2name)
is_element_regex<-function(term,regs){return(any(sapply(regs,grepl,x=term,ignore.case=T)))}
bblockers_in_data = medcode2name[sapply(medcode2name,is_element_regex,regs=bblockers)]
bblockers_in_data = name2medcode[bblockers_in_data]
ccblockers_in_data = medcode2name[sapply(medcode2name,is_element_regex,regs=ccblockers)]
ccblockers_in_data = name2medcode[ccblockers_in_data]
is_treatment_in_set<-function(x,s){
  x = x[!is.na(x)]
  return(any(is.element(x,set=s)))
}
sample2bblocker = apply(medication_data,1,is_treatment_in_set,s=bblockers_in_data)
table(sample2bblocker)
sample2ccblocker = apply(medication_data,1,is_treatment_in_set,s=ccblockers_in_data)
table(sample2ccblocker)
# NO DEATH DATA AT THIS POINT
#death_cols = colnames(pheno_data)[grepl(pattern="death",colnames(pheno_data),ignore.case = T)]
rm(pheno_data);gc()

# Get ICD data
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
icd_data_code_classes[["bradycardia"]] = "R001"
icd_data_code_classes[["tachycardia"]] = "R000"
icd_data_code_columns = lapply(icd_data_code_classes,intersect,y=colnames(icd_data))
icd_data_code_columns = icd_data_code_columns[sapply(icd_data_code_columns,length)>0]
get_icd_data<-function(x,y){
  if(length(x)==1){v=y[,x];names(v)=rownames(y);return(v)}
  return(y[,x])
}
icd_outcome_matrices = lapply(icd_data_code_columns,get_icd_data,y=icd_data)
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

# Define the data for the observational analysis
subjects_with_data = intersect(rownames(simple_covs),euro_ids) # European that have covariates
subjects_with_data = intersect(subjects_with_data,rownames(geno_data)) # Have genotypes
subjects_with_data = intersect(subjects_with_data,rownames(geno_data)) # In icd data
reduce_to_subjs<-function(x,y){
  if(is.null(dim(x))){return(x[y])}
  return(x[y,])
}
covariate_matrix = covariate_matrix[subjects_with_data,]
icd_data = icd_data[subjects_with_data,]
icd_outcome_matrices = lapply(icd_outcome_matrices,reduce_to_subjs,y=subjects_with_data)
geno_data = geno_data[subjects_with_data,]
medication_data = medication_data[subjects_with_data,]
additional_scores = lapply(additional_scores,reduce_to_subjs,y=subjects_with_data)
simple_covs = simple_covs[subjects_with_data,]
sample2bblocker = sample2bblocker[subjects_with_data]
sample2ccblocker = sample2ccblocker[subjects_with_data]
gc()

# load MR results
load("/Users/David/Desktop/ukbb/gwas_interpretation/mr/results_nov_1_2017_slim.RData")

##################################################
##################################################
##################################################
# Save the current workspace - should make it easier to reload and analyze
save(list = ls(all.names = TRUE), file = "/Users/David/Desktop/ukbb/gwas_interpretation/observational_analysis_workspace.RData", envir = .GlobalEnv)
##################################################
##################################################
##################################################
load( "/Users/David/Desktop/ukbb/gwas_interpretation/observational_analysis_workspace.RData")
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
outcome_vs = lapply(icd_outcome_matrices,function(x,...)rowSums(as.matrix(x),...),na.rm=T)
outcome_vs[1:4] = lapply(outcome_vs[1:4],transform_outcome_to_binary)

expo_vs = list()
expo_vs[["ExerciseHR"]] = expo[,2]
expo_vs[["Recovery"]] = expo[,1]
expo_vs = lapply(expo_vs,function(x)x[!is.na(x)])
expo_vs = lapply(expo_vs,gaus_norm)

# Create a data frame with data for subsequent analyses
heart_disease_v = outcome_vs$heart_disease
is_healthy = is.element(names(heart_disease_v),set=healthy_subjs)
is_healthy = as.numeric(is_healthy);names(is_healthy)=names(heart_disease_v)
AGE = simple_covs[,2];names(AGE) = rownames(simple_covs)
SEX = simple_covs[,1];names(SEX) = rownames(simple_covs)
colnames(covariate_matrix)[grepl(pattern = "chol",colnames(covariate_matrix),ignore.case = T)]
RBC = covariate_matrix[,"Red blood cell (erythrocyte) count"];names(RBC)=rownames(covariate_matrix)
HEMOC = covariate_matrix[,"Haemoglobin concentration"];names(HEMOC)=rownames(covariate_matrix)
HEMOMCH = covariate_matrix[,"Mean corpuscular haemoglobin concentration"];names(HEMOMCH)=rownames(covariate_matrix)
FAT = covariate_matrix[,"Body fat percentage"];names(FAT)=rownames(covariate_matrix)
HEIGHT = covariate_matrix[,"Standing height"];names(HEIGHT) = rownames(covariate_matrix)
INCOME = covariate_matrix[,"Average total household income before tax"];names(INCOME) = rownames(covariate_matrix)
MEDI = medication_data[,"Medication for cholesterol, blood pressure, diabetes, or take exogenous hormones.0.0"]
SMOKING = covariate_matrix[,"Current tobacco smoking"];names(SMOKING)=rownames(covariate_matrix)
BRADY = outcome_vs$bradycardia
TACHY = outcome_vs$tachycardia
MI_Disease = as.numeric(rowSums(icd_data[,c("I219","I210","I211","I214")],na.rm=T)>0);names(MI_Disease)=rownames(icd_data)
BBLOCKERS = sample2bblocker
CCBLOCKERS = sample2ccblocker
HeartDISEASE = heart_disease_v
NUM_DISEASES = outcome_vs$disease_merge
HYPERTENSION = outcome_vs$hypertension
DIABETES = outcome_vs$diabetes
CANCER = outcome_vs$cancer
DBP = covariate_matrix[,"Diastolic blood pressure, automated reading"];names(DBP) = rownames(covariate_matrix)
SBP = covariate_matrix[,"Systolic blood pressure, automated reading"];names(SBP) = rownames(covariate_matrix)
inds = intersect(subjects_with_data,names(expo_vs$ExerciseHR))
EHR = expo_vs$ExerciseHR[inds] # Specific to our tests
inds = intersect(subjects_with_data,names(expo_vs$Recovery))
REC = expo_vs$Recovery[inds] # Specific to our tests
RHR = additional_scores$`Pulse rate`[subjects_with_data]
POLLUTION = covariate_matrix[,"Nitrogen dioxide air pollution; 2007"];names(POLLUTION) = rownames(covariate_matrix)
IN_EXERCISE_COHORT = is.element(subjects_with_data,set=union(names(EHR),names(REC)));names(IN_EXERCISE_COHORT)=subjects_with_data
AFIB = icd_data[,"I48"];names(AFIB)=rownames(icd_data)
MAXWD = fitness_scores_matrix[intersect(names(EHR),names(REC)),"Max achieved WD"]
hf_cols = colnames(icd_data)[sapply(colnames(icd_data),has_any_regex,regexs = c("I42","I50"))]
HF = as.numeric(rowSums(icd_data[,hf_cols],na.rm=T)>0);names(MI_Disease)=rownames(icd_data)
# Add "Trunk fat-free mass"
TRUNK_MASS = covariate_matrix[,"Trunk fat-free mass"];names(TRUNK_MASS)=rownames(covariate_matrix)
# Add "Leg fat-free mass (right)"
LEG_FATFREE_MASS = covariate_matrix[,"Leg fat-free mass (right)"];names(LEG_FATFREE_MASS)=rownames(covariate_matrix)
# Add "Leg fat percentage (right)"
LEG_FAT_PER = covariate_matrix[,"Leg fat percentage (right)"];names(LEG_FAT_PER)=rownames(covariate_matrix)
pairwise_cor(FAT,TRUNK_MASS)
pairwise_cor(FAT,LEG_FAT_PER)
pairwise_cor(LEG_FATFREE_MASS,TRUNK_MASS)
# Add "Lymphocyte count"
LYMPHC = covariate_matrix[,"Lymphocyte count"];names(LYMPHC)=rownames(covariate_matrix)
# Add "Basal metabolic rate"
METABOLIC_RATE = covariate_matrix[,"Basal metabolic rate"];names(METABOLIC_RATE)=rownames(covariate_matrix)
# Add "Body mass index (BMI)"
BMI = covariate_matrix[,"Body mass index (BMI)"]; names(BMI)=rownames(covariate_matrix)
pairwise_cor(BMI,FAT)
pairwise_cor(BMI,TRUNK_MASS)
pairwise_cor(BMI,METABOLIC_RATE)
# Add "Waist circumference"
WAISTC = covariate_matrix[,"Waist circumference"];names(WAISTC)=rownames(covariate_matrix)
pairwise_cor(BMI,WAISTC)
# Add "Alcohol drinker status" or "Alcohol intake frequency"
ALCOHOL1 = 7-as.numeric(covariate_matrix[,"Alcohol intake frequency"]);names(ALCOHOL1)=rownames(covariate_matrix)
# Add "Frequency of depressed mood in last 2 weeks"  or "Frequency of unenthusiasm / disinterest in last 2 weeks" 
#     or "Frequency of tenseness / restlessness in last 2 weeks"
MOOD_MAT = covariate_matrix[,c("Frequency of depressed mood in last 2 weeks","Frequency of unenthusiasm / disinterest in last 2 weeks",
                               "Frequency of tenseness / restlessness in last 2 weeks")]
MOOD_MAT = sapply(MOOD_MAT,as.numeric) # This matrix is highly correlated, summarize in a single score
MOOD = apply(MOOD_MAT,1,sum,na.rm=T);names(MOOD) = rownames(covariate_matrix)
# Look at "Coffee intake"
# COFFEE = covariate_matrix[,"Coffee intake"]
# COFFEE[COFFEE==-10]=0.5 # Change less than 1 to 0.5 (did not make sense...)
# Look at "Sleep duration"
SLEEP = covariate_matrix[,"Sleep duration"];names(SLEEP) = rownames(covariate_matrix)
# Look at "Time spent using computer"
# Add Spirometry data
SPIRO_FEV1 = apply(
  additional_scores$Spirometry[,grepl(colnames(additional_scores$Spirometry),pattern="Forced expiratory volume in 1-second \\(FEV1\\)")],
  1,mean,na.rm=T)
names(SPIRO_FEV1) = rownames(additional_scores$Spirometry)
SPIRO_FVC = apply(
  additional_scores$Spirometry[,grepl(colnames(additional_scores$Spirometry),pattern="Forced vital capacity \\(FVC\\)")],
  1,mean,na.rm=T)
names(SPIRO_FVC) = rownames(additional_scores$Spirometry)
SPIRO_PEF = apply(
  additional_scores$Spirometry[,grepl(colnames(additional_scores$Spirometry),pattern="Peak expiratory flow \\(PEF\\)")],
  1,mean,na.rm=T)
names(SPIRO_PEF) = rownames(additional_scores$Spirometry)
pairwise_cor(SPIRO_FVC,SPIRO_PEF)
pairwise_cor(SPIRO_FVC,SPIRO_FEV1)
pairwise_cor(SPIRO_FVC,DATA$RBC)
# Hand grip strength
HAND_GRIP = apply(
  additional_scores$`hand grip`[,grepl(colnames(additional_scores$`hand grip`),pattern="right")],
  1,mean,na.rm=T)

# Get the genetic score
# Get loadings using multivariate analysis, correct for AGE, SEX, BBLOCKERS, PCs
snp_lists = lapply(snp_lists,intersect,y=colnames(geno_data))
PCs = covariate_matrix[,grepl("PC",colnames(covariate_matrix))]
inds = names(expo_vs$ExerciseHR)
inds = intersect(inds,subjects_with_data)

# Inspect a specific genetic weighted score: heart disease has inverted results
nn = names(mr_analysis_res)[1]
bx = mr_analysis_res[[nn]]$multivar_input@betaX
mm = as.matrix(geno_data[,names(bx)])
mm[is.na(mm)] = 0
gws1 = mm %*% bx
names(gws1) = rownames(geno_data)
rm(mm);gc()
GWS=gws1

DATA = data.frame(AGE,SEX,FAT,HEIGHT,INCOME,BRADY,
                  TACHY,BBLOCKERS,CCBLOCKERS,HeartDISEASE,MI_Disease,
                  DIABETES,SBP,DBP,NUM_DISEASES,HYPERTENSION,
                  POLLUTION,IN_EXERCISE_COHORT,GWS,RHR,AFIB,HF,
                  CANCER,SMOKING,HEMOC,HEMOMCH,RBC,
                  TRUNK_MASS,LYMPHC,METABOLIC_RATE,BMI,
                  MOOD,SLEEP,SPIRO_FEV1,SPIRO_FVC,HAND_GRIP)
# Some plots and comparisons
boxplot(GWS~as.numeric(IN_EXERCISE_COHORT),data=DATA,las=2)
boxplot(AGE~as.numeric(IN_EXERCISE_COHORT),data=DATA,las=2)
summary(glm(IN_EXERCISE_COHORT~.,data=DATA,family = binomial(link='logit')))
summary(lm(SBP~HEIGHT+GWS+SMOKING+RBC,data=DATA))

# Split snps and regions using gwascat data
region_to_traits = list()
region_to_traits[["ExerciseHR"]] = unique(cbind(gwascat_data$ExerciseHR$GenomicLocus,gwascat_data$ExerciseHR$IndSigSNP,
                                                gwascat_data$ExerciseHR$Trait))
unmapped_regions_ehr = setdiff(region_data$ExerciseHR$GenomicLocus,gwascat_data$ExerciseHR$GenomicLocus)
region_to_traits[["Recovery"]] = unique(cbind(gwascat_data$Recovery$GenomicLocus,gwascat_data$Recovery$IndSigSNP,
                                                gwascat_data$Recovery$Trait))
lapply(region_to_traits,function(x)sort(table(x[,3]),decreasing=T)[1:10])
unmapped_regions_rec = setdiff(region_data$Recovery$GenomicLocus,gwascat_data$Recovery$GenomicLocus)
geno_data_snp_info = c()
for(snp in colnames(geno_data)){
  regions = lapply(region_data,function(x,y)x[x$rsID==y,1],y=snp)
  traits = c()
  for(i in 1:length(regions)){
    r = regions[[i]]
    if(length(r)==0){next}
    m = region_to_traits[[i]]
    traits = c(traits,m[m[,1]==r,3])
  }
  geno_data_snp_info[[snp]] = traits
}

# Try to learn a causal model from the data
discDATA = lapply(DATA,disc_data_using_cut,cuts=10)
discDATA = data.frame(discDATA)
rownames(discDATA)=rownames(DATA)
save(DATA,discDATA,covariate_matrix,PCs,geno_data,EHR,REC,expo_vs,snp_lists,MAXWD,MEDI,
     file = "/Users/David/Desktop/ukbb/gwas_interpretation/observational_analysis_data.RData", envir = .GlobalEnv)
# # Test 
# dd = DATA[,c("HEIGHT","SEX","RHR","DBP")]
# dd = dd[!apply(is.na(dd),1,any),]
# cor.test(lm(HEIGHT~SEX+RHR,dd)$residuals,lm(DBP~SEX+RHR,dd)$residuals)$p.value
# run_ci_test_one_is_numeric(1,4,2:3,dd)
# # Check cond indep between rhr and gws scores
# for(nn1 in colnames(discDATA)){
#   if(nn1=="RHR" || grepl(pattern="GWS",nn1)){next}
#   for(nn2 in colnames(discDATA)){
#     if(nn2=="RHR" || nn1==nn2 || grepl(pattern="GWS",nn2)){next}
#     p1 = run_discrete_ci_test("RHR","GWS",c(nn1,nn2),discDATA)
#     p2 = run_discrete_ci_test("RHR","GWSHeight",c(nn1,nn2),discDATA)
#     p3 = run_discrete_ci_test("RHR","GWSnonheight",c(nn1,nn2),discDATA)
#     print(c(c(nn1,nn2),p1,p2,p3))
#   }
# }

# Attempt 4: analyze the SNPs in the exercise cohort
# An alternative analysis of the detected SNPs vs. the selected phenotypes of all cohorts
TRAITS = c("AGE","SEX","FAT","HEIGHT","RHR","DBP","SBP","BRADY","RBC","HEMOMCH","BBLOCKERS","CCBLOCKERS",
           "LYMPHC","BMI","MOOD","SLEEP","TRUNK_MASS","METABOLIC_RATE","SPIRO_FEV1","SPIRO_FVC",
           "HAND_GRIP","DIABETES","HeartDISEASE")
#setdiff(TRAITS,colnames(DATA))
length(TRAITS)
subjs = intersect(names(EHR),names(REC))
subjs = intersect(subjs,rownames(DATA))
snp_analysis_data = data.frame(discDATA[subjs,TRAITS],
                               EHR=disc_data_using_cut(EHR[subjs],cuts=10),
                               REC=disc_data_using_cut(REC[subjs],cuts=10))
snp_analysis_contindata = data.frame(DATA[subjs,TRAITS],
                               EHR=EHR[subjs],REC=REC[subjs])

summary(lm(EHR~.,data=snp_analysis_contindata))
qqnorm(scale(lm(EHR~.,data=snp_analysis_contindata)$residuals))
abline(0,1)
summary(lm(REC~.,data=snp_analysis_contindata))
qqnorm(scale(lm(REC~.,data=snp_analysis_contindata)$residuals))
abline(0,1)

TRAITS = c(TRAITS,"EHR","REC")
snp_analysis_data = snp_analysis_data[!apply(is.na(snp_analysis_data),1,any),]
snp_analysis_contindata = snp_analysis_contindata[!apply(is.na(snp_analysis_contindata),1,any),]
gc()
snp_vs_trait_ci_network = c()
sepsets = list()
for (i in 1:ncol(geno_data)){
  currsnp = colnames(geno_data)[i]
  sepsets[[currsnp]]=list()
  snpv = factor(geno_data[rownames(snp_analysis_data),currsnp])
  inds = !is.na(snpv)
  gc()
  currdata = list(DATA=data.frame(snp=snpv[inds],snp_analysis_contindata[inds,]),
                  discDATA=data.frame(snp=factor(snpv[inds]),snp_analysis_data[inds,]))
  curredges = c()
  for (target in 4:26){
    has_edge1 = analyze_source_vs_target_node(1,target,c(2:3),currdata,run_ci_test,depth = 2,p0 = 0.05,pthr = 0.05)
    has_edge2 =  analyze_source_vs_target_node(1,target,NULL,currdata,run_ci_test,depth = 2,p0 = 0.05,pthr = 0.05)
    sepset = list(colnames(currdata$DATA)[has_edge1$sepset],colnames(currdata$DATA)[has_edge2$sepset])
    has_edge=has_edge1$res & has_edge2$res
    print(c(currsnp,colnames(currdata$DATA)[target],has_edge))
    if(any(sapply(sepset,length)>0)){print(c(unlist(sepset)))}
    curredges = c(curredges,has_edge)
    sepsets[[currsnp]][[colnames(currdata$DATA)[target]]]=sepset
  }
  names(curredges) = colnames(currdata$DATA)[4:26]
  snp_vs_trait_ci_network = rbind(snp_vs_trait_ci_network,curredges)
  rownames(snp_vs_trait_ci_network)[nrow(snp_vs_trait_ci_network)]=currsnp
  print(paste("DONE with snp number ",i, "id:", currsnp))
}
table(rowSums(snp_vs_trait_ci_network))
sort(rowSums(snp_vs_trait_ci_network))
colSums(snp_vs_trait_ci_network)
apply(snp_vs_trait_ci_network,2,function(x)names(which(x)))
snp2traits = apply(snp_vs_trait_ci_network,1,function(x)names(which(x)))
table(rowSums(snp_vs_trait_ci_network)==0)
table(rowSums(snp_vs_trait_ci_network))
exercise_genetic_scores=c()
bx_ehr = mr_analysis_res[[1]]$multivar_input@betaX
bx_rec = mr_analysis_res[[13]]$multivar_input@betaX
for(j in 1:ncol(snp_vs_trait_ci_network)){
  curr_snps = rownames(snp_vs_trait_ci_network)[snp_vs_trait_ci_network[,j]]
  if(length(curr_snps)< 5){next}
  bx1 = bx_ehr[curr_snps]
  bx1[is.na(bx1)]=0
  bx2 = bx_rec[curr_snps]
  bx2[is.na(bx2)] = 0
  currbx = bx1+bx2;names(currbx)=curr_snps
  mm = as.matrix(geno_data[,names(currbx)])
  mm[is.na(mm)] = 0
  gws1 = mm %*% currbx
  exercise_genetic_scores = cbind(exercise_genetic_scores,gws1)
  colnames(exercise_genetic_scores)[ncol(exercise_genetic_scores)] = colnames(snp_vs_trait_ci_network)[j]
}
cor(exercise_genetic_scores[rownames(DATA),],method="spearman")

ENV = c("INCOME","POLLUTION","SMOKING")

save(DATA,discDATA,covariate_matrix,PCs,gws1,gws2,geno_data,EHR,REC,expo_vs,snp_lists,MAXWD,MEDI,
     TRAITS,ENV,snp_analysis_data,snp_analysis_contindata,exercise_genetic_scores,snp_vs_trait_ci_network,
     file = "/Users/David/Desktop/ukbb/gwas_interpretation/observational_analysis_data.RData", envir = .GlobalEnv)

# Genetic scores: compute correlations with diseases and EHR
EXSCORES = exercise_genetic_scores[rownames(DATA),]
colnames(EXSCORES) = paste("GS",colnames(EXSCORES),sep="_")
DATA = data.frame(DATA,EXSCORES)
summary(glm(HeartDISEASE~GS_RHR+SEX+AGE,data=DATA,family = binomial(link='logit')))
summary(glm(HeartDISEASE~GS_REC+GS_RHR+SEX+AGE,data=DATA,family = binomial(link='logit')))
summary(lm(HEIGHT~GS_HEIGHT+SEX+AGE,data=DATA))

# # MR analysis of RHR 
# ehr_snps = rownames(snp_vs_trait_ci_network)[snp_vs_trait_ci_network[,"EHR"]]
# ehr_snp_data = geno_data[,ehr_snps]
# ehr_snp_data[is.na(ehr_snp_data)]=0
# RHR = DATA$RHR;names(RHR)=rownames(DATA)
# ehr_rhr_mr =  run_mr_analysis(EHR,RHR,ehr_snp_data,covs = DATA[,c("SEX","AGE")])
# mr_plot(ehr_rhr_mr$multivar_input)
# rec_snps = rownames(snp_vs_trait_ci_network)[snp_vs_trait_ci_network[,"REC"]]
# rec_snp_data = geno_data[,rec_snps]
# rec_snp_data[is.na(rec_snp_data)]=0
# rec_rhr_mr =  run_mr_analysis(REC,RHR,rec_snp_data,covs = DATA[,c("SEX","AGE")])
# mr_plot(rec_rhr_mr$multivar_input)
# rec_rhr_mr$multivar_all

DISEASES = c("HeartDISEASE","HYPERTENSION","DIABETES","CANCER")
for(dis in DISEASES){
  print(summary(glm(as.formula(paste(dis,"~","AGE+SEX+GWS")),data=DATA,family = binomial(link='logit'))))
}
glm_objects=list()
for(dis in DISEASES){
  print(dis)
  glm0 = glm(as.formula(paste(dis,"~","AGE+SEX")),data=DATA,family = binomial(link='logit'))
  single_glms = list()
  for(gs in colnames(EXSCORES)){
    single_glms[[gs]] = glm(as.formula(paste(dis,"~",gs,"+AGE+SEX")),data=DATA,family = binomial(link='logit'))
  }
  all_glm = glm(as.formula(paste(dis,"~",paste(colnames(EXSCORES),collapse="+"),"+AGE+SEX")),data=DATA,family = binomial(link='logit'))
  glm_objects[[dis]]=list()
  glm_objects[[dis]][["0"]] = glm0
  glm_objects[[dis]][["singles"]] = single_glms
  glm_objects[[dis]][["all"]] = all_glm
}
# Interpret the results: single scores
anova_results = list()
for(dis in DISEASES){
  anova_results[["ds"]] = list()
  for(gs in colnames(EXSCORES)){
    av = anova(glm_objects[[dis]][["0"]],glm_objects[[dis]][["singles"]][[gs]],test="LRT")
    anova_results[[dis]][[gs]]=av
  }
}
# Get the best single score and compare to using all
model_selection_ps = c()
all_single_ps = c()
for(dis in DISEASES){
  ps = sapply(anova_results[[dis]],function(x)x$`Pr(>Chi)`[2])
  all_single_ps = rbind(all_single_ps,ps)
  print(dis)
  print(p.adjust(ps,method='fdr'))
  min_i = which(ps==min(ps))[1]
  gs = names(ps)[min_i]
  all_p = anova(glm_objects[[dis]][["singles"]][[gs]],glm_objects[[dis]][["all"]],test="LRT")
  model_selection_ps = rbind(model_selection_ps,
                             c(dis,gs,ps[min_i],all_p$`Pr(>Chi)`[2]))
}
rownames(all_single_ps) = DISEASES
allps = c(all_single_ps)
allqs = p.adjust(allps,method='fdr')
thr = max(allps[allqs<0.1])
gs2disease=apply(all_single_ps,2,function(x,y)names(which(x<=y)),y=thr)
names(gs2disease) = sapply(names(gs2disease),gsub,pattern="GS_",replace="")

# Print the clustering analysis results for display items and interpretation
# Fuma
snp_locs = read.delim("/Users/David/Desktop/ukbb/gwas_interpretation/id.txt",stringsAsFactors = F,sep=" ")
rownames(snp_locs) = snp_locs[,3]
for(j in 1:ncol(snp_vs_trait_ci_network)){
  curr_snps = rownames(snp_vs_trait_ci_network)[snp_vs_trait_ci_network[,j]]
  if(length(curr_snps)< 5){next}
  bx1 = bx_ehr[curr_snps]
  bx1[is.na(bx1)]=0
  bx2 = bx_rec[curr_snps]
  bx2[is.na(bx2)] = 0
  print(c(colnames(snp_vs_trait_ci_network)[j],sum(abs(bx1)),sum(abs(bx2))))
  m = cbind(curr_snps,rep(1e-10,length(curr_snps)))
  colnames(m) = c("rsID","P-value")
  fname = paste("/Users/David/Desktop/ukbb/gwas_interpretation/fuma_in_files/interpretation_",
                colnames(snp_vs_trait_ci_network)[j],".txt",sep='')
  write.table(m,file=fname,row.names = F,col.names = T,quote = F,sep="\t")
}
# SNP supp table
bx_ehr = mr_analysis_res[[1]]$multivar_input@betaX
bx_rec = mr_analysis_res[[13]]$multivar_input@betaX
# Trait network for the analysis
edges = c()
for(j in 1:ncol(snp_vs_trait_ci_network)){
  curr_snps = rownames(snp_vs_trait_ci_network)[snp_vs_trait_ci_network[,j]]
  if(length(curr_snps)< 5){next}
  curr_name = colnames(snp_vs_trait_ci_network)[j]
  curr_ds = gs2disease[[curr_name]]
  if(curr_name=="EHR"){
    curr_name = "ExerciseHR"
  }
  if(curr_name=="Rec"){
    currname = "Recovery"
  }
  curr_name = gsub(curr_name,pattern="SPIRO_",replace="")
  curr_name = paste(curr_name,'(',length(curr_snps),')',sep='')
  print(curr_name)
  bx1 = bx_ehr[curr_snps]
  bx1[is.na(bx1)]=0
  bx2 = bx_rec[curr_snps]
  bx2[is.na(bx2)] = 0
  We = sum(abs(bx1))
  Wr = sum(abs(bx2))
  edges = rbind(edges,c(curr_name,"ExerciseHR",We))
  edges = rbind(edges,c(curr_name,"Recovery",Wr))
  for(d in curr_ds){
    edges = rbind(edges,c(curr_name,d,1))
  }
}
edges = edges[as.numeric(edges[,3])>0,]
write.table(edges,file="/Users/David/Desktop/ukbb/gwas_interpretation/causal_inference/traits_network_v2.txt",
            row.names = F,col.names = F,quote=F,sep="\t")

####### Load mortality data and run MRs#########
mortality_data = read.delim('/Users/David/Desktop/ukbb/gwas_interpretation/ukb11180.tab')
dim(mortality_data)
table(mortality_data[,"f.40001.0.0"])
# Exclude mortalities with ICD starting with R,S,T,V,Y
primary_mortality_data = mortality_data[,c("f.40001.0.0","f.40001.1.0","f.40001.2.0")]
primary_mortality_data = as.matrix(primary_mortality_data)
primary_mortality_data[is.na(primary_mortality_data)]="0"
all_death_cases = apply(primary_mortality_data,1,function(x)sum(x!="0")>0) 
table(all_death_cases)
excluded_mortalities = apply(primary_mortality_data,1,function(x)any(sapply(x,grepl,pattern="^(R|S|T|V|Y)"))) 
cancer_mortalities = apply(primary_mortality_data,1,function(x)any(sapply(x,grepl,pattern="^(C|D)"))) 
table(excluded_mortalities)
table(cancer_mortalities)
MORTALITY = as.numeric(all_death_cases & !excluded_mortalities & !cancer_mortalities)
names(MORTALITY) = mortality_data[,1]
table(MORTALITY[names(EHR)])
ehr_mort_glm = glm(factor(MORTALITY)~.,data=data.frame(MORTALITY=MORTALITY[names(EHR)],
                                                       SEX=DATA[names(EHR),"SEX"],EHR=EHR,
                                                       AGE=DATA[names(EHR),"AGE"]),family = binomial(link='logit'))
summary(ehr_mort_glm)
AGE_OF_DEATH = apply(mortality_data[,c("f.40007.0.0","f.40007.1.0","f.40007.2.0")],1,mean,na.rm=T)
table(is.nan(AGE_OF_DEATH))
names(AGE_OF_DEATH) = names(MORTALITY)
TIME_TO_DEATH = AGE_OF_DEATH[rownames(DATA)] - DATA$AGE
MORTALITY = MORTALITY[rownames(DATA)]
table(TIME_TO_DEATH>5)
table(TIME_TO_DEATH,MORTALITY)
MORTALITY[TIME_TO_DEATH>5] = 0
ehr_mort_glm = glm(factor(MORTALITY)~.,data=data.frame(MORTALITY=MORTALITY[names(EHR)],
      SEX=DATA[names(EHR),"SEX"],EHR=EHR,
      AGE=DATA[names(EHR),"AGE"]),family = binomial(link='logit'))
summary(ehr_mort_glm)
###########################################################


# MR for each one of our scores
library(MendelianRandomization)
mortality_mr_res = list()
OV = MORTALITY[!is.na(names(MORTALITY))]
for (gs in colnames(EXSCORES)){
  trait = gsub(gs,pattern="GS_",replace="")
  trait = "EHR"
  curr_snps = rownames(snp_vs_trait_ci_network)[snp_vs_trait_ci_network[,trait]]
  curr_d = geno_data[rownames(DATA),curr_snps]
  table(colSums(is.na(curr_d)))
  if(trait=="EHR"){
    TR = EHR
  }
  if(trait=="REC"){
    TR = REC
  }
  if(trait!="REC" && trait!="EHR"){
    TR = DATA[,trait]
    names(TR) = rownames(DATA)
  }
  covs = DATA[names(TR),c("SEX","AGE")]
  mortality_mr_res[[gs]] = run_mr_analysis(expo_v = TR,outcome_v = OV,snps = curr_d,covs = covs)
  mr_plot(mortality_mr_res[[gs]]$multivar_input)
  print(mortality_mr_res[[gs]]$multivar_egger)
}

# MR of the original scores
summary(lm(MORTALITY~.,data=data.frame(MORTALITY=OV,DATA[names(OV),c("GWSHeight","SEX","AGE")])))
summary(lm(MORTALITY~.,data=data.frame(MORTALITY=OV,DATA[names(OV),c("GWSnonheight","SEX","AGE")])))
summary(lm(MORTALITY~.,data=data.frame(MORTALITY=OV,DATA[names(OV),c("GWS","SEX","AGE")])))
summary(lm(MORTALITY~.,data=data.frame(MORTALITY=OV,DATA[names(OV),c("GS_AGE","SEX","AGE")])))
summary(lm(MORTALITY~.,data=data.frame(MORTALITY=OV,DATA[names(OV),c("GS_EHR","SEX","AGE")])))
summary(lm(MORTALITY~.,data=data.frame(MORTALITY=OV,DATA[names(OV),c("GS_RHR","SEX","AGE")])))
summary(lm(MORTALITY~.,data=data.frame(MORTALITY=OV,DATA[names(OV),c("GS_HEIGHT","SEX","AGE")])))
summary(lm(MORTALITY~.,data=data.frame(MORTALITY=OV,DATA[names(OV),])))
lm1 = glm(factor(MORTALITY)~.,data=data.frame(MORTALITY=OV,DATA[names(OV),c("GS_REC","GS_EHR","GS_RHR","GS_HEIGHT","GS_AGE","SEX","AGE")]),
            family = binomial(link='logit'))
lm2 = glm(factor(MORTALITY)~.,data=data.frame(MORTALITY=OV,DATA[names(OV),c("SEX","AGE")]),
          family = binomial(link='logit'))
anova(lm1,lm2,test="LRT")

summary(glm(factor(SMOKING)~.,data=data.frame(DATA[,c("GS_REC","GS_EHR","GS_RHR","GS_HEIGHT","GS_AGE","SEX","AGE","SMOKING")]),
          family = binomial(link='logit')))

summary(lm(GS_EHR~RHR+SEX+AGE,data=DATA[names(EHR),]))
summary(lm(GS_EHR~RHR+SEX+AGE+EHR,data=data.frame(EHR=EHR,DATA[names(EHR),])))
summary(lm(GS_RHR~RHR+SEX+AGE+EHR,data=data.frame(EHR=EHR,DATA[names(EHR),])))

# # tests
# snpv = geno_data[rownames(snp_analysis_data),1]
# inds = !is.na(snpv)
# dd = list(DATA=data.frame(snp=snpv[inds],snp_analysis_contindata[inds,]),
#           discDATA=data.frame(snp=factor(snpv[inds]),snp_analysis_data[inds,]))
# run_ci_test(1,"EHR",c(2:3,12),data=dd)
# summary(lm(EHR~snp+AGE+SEX+PC1+PC2,data=dd$DATA))
# run_discrete_ci_test(1,"EHR",c(2:3,12),data=dd$discDATA,test='mi-cg')
# dd2 = apply(dd$DATA,2,disc_data_using_cut,cuts=5)
# dd2 = list()
# for(j in 1:ncol(dd$DATA)){
#   dd2[[colnames(dd$DATA)[j]]] = disc_data_using_cut(dd$DATA[,j],cuts=3)
#   print(j)
# }
# dd2 = as.data.frame(dd2)
# colnames(dd2) = colnames(dd$DATA);rownames(dd2) = rownames(dd$DATA)
# apply(dd2,2,function(x)length(unique(x)))
# run_discrete_ci_test(1,"EHR",c(2:3,12),data=dd2,test='mi')

# Try causal inference packages
load("/Users/David/Desktop/ukbb/gwas_interpretation/observational_analysis_data.RData")
DATA$BBLOCKERS = factor(DATA$BBLOCKERS)
DATA$CCBLOCKERS = factor(DATA$CCBLOCKERS)
DATA$SEX = factor(DATA$SEX)
DATA$AFIB = factor(DATA$AFIB)
DATA$MI_Disease = factor(DATA$MI_Disease)
DATA$TACHY = factor(DATA$TACHY)
DATA$BRADY = factor(DATA$BRADY)
DATA$HeartDISEASE = factor(DATA$HeartDISEASE)
DATA$DIABETES = factor(DATA$DIABETES)
DATA$HYPERTENSION = factor(DATA$HYPERTENSION)
DATA$AGE = as.numeric(DATA$AGE)
DATA$SBP = as.numeric(DATA$SBP)
DATA$DBP = as.numeric(DATA$DBP)
names(which(sapply(DATA,is.numeric)))
names(which(sapply(DATA,is.integer)))
inds = c("AGE","SEX","RHR","GWSHeight","GWSnonheight",
         "DIABETES","AFIB","HEIGHT","FAT","SBP",
         "DBP","MI_Disease","INCOME","SMOKING","BRADY",
         "BBLOCKERS","CCBLOCKERS","POLLUTION")
table(sapply(DATA[,inds],is.numeric))
# USE pcalg to get a skeleton
causal_skeleton_nmax5_e20 = pcalg::skeleton(list(DATA=DATA[,inds],discDATA=discDATA[,inds]), 
                    run_ci_test, alpha = 1e-20, verbose=TRUE,labels=inds,m.max=5)
causal_skeleton_nmax5_e50 = pcalg::skeleton(list(DATA=DATA[,inds],discDATA=discDATA[,inds]),
                    run_ci_test, alpha = 1e-50, verbose=TRUE,labels=inds,m.max=5)
causal_skeleton_nmax5_e10 = pcalg::skeleton(list(DATA=DATA[,inds],discDATA=discDATA[,inds]), 
                    run_ci_test, alpha = 1e-10, verbose=TRUE,labels=inds,m.max=5)
causal_skeleton_nmax5 = pcalg::skeleton(list(DATA=DATA[,inds],discDATA=discDATA[,inds]), 
                    run_ci_test, alpha = 1e-5, verbose=TRUE,labels=inds,m.max=5)

save(causal_skeleton_nmax5,causal_skeleton_nmax5_e10,causal_skeleton_nmax5_e20,causal_skeleton_nmax5_e50,file=
       "/Users/David/Desktop/ukbb/gwas_interpretation/causal_inference/skeletal_graphs_nov13_2017.RData")
plot(causal_skeleton_nmax5)
plot(causal_skeleton_nmax5_e10)
dseps = causal_skeleton_nmax5_e10@sepset

# AFIB vs. Height
i1 = 7; i2= 8
inds[c(i1,i2)]
sepset = union(dseps[[i1]][[i2]],dseps[[i2]][[i1]])
inds[sepset]
run_ci_test(i1,i2,sepset,list(DATA=DATA[,inds],discDATA=discDATA[,inds]))
run_discrete_ci_test(i1,i2,sepset,discDATA[,inds])

gaps = 1-as(causal_skeleton_nmax5@graph,"matrix")
table(gaps)

# rinds = !apply(is.na(DATA[,c("AGE","RHR","HEIGHT")]),1,any)
# ci.test("AGE","HEIGHT","RHR",data=DATA[rinds,c("AGE","RHR","HEIGHT")],test="cor")$p.value
# run_ci_test_one_is_numeric("AGE","HEIGHT","RHR",data=DATA[rinds,])
causal_model = rfci(list(DATA=DATA[,inds],discDATA=discDATA[,inds]), run_ci_test, alpha = 1e-5, verbose=TRUE,fixedGaps = gaps,
                   labels=inds,m.max = 0,conservative = F)
plot(causal_model)
show(causal_model)
slot(causal_model,"allPdsep")

# # Estimate triplets based on the CIT package
# library(MendelianRandomization)
# triplet = inds[c(4,8,7)]
# covnames = c("SEX","AGE","RHR")
# rinds = !apply(is.na(DATA[,union(covnames,triplet)]),1,any)
# covs = DATA[rinds,covnames]
# Lv = DATA[rinds,triplet[1]]
# Gv = DATA[rinds,triplet[2]]
# Tv = as.numeric(as.character(DATA[rinds,triplet[3]]))
# Rg = lm(Gv~.,data=data.frame(Gv=Gv,covs))$residuals
# Rt = lm(Tv~.,data=data.frame(Tv=Tv,covs))$residuals
# bx = summary(lm(Rg~Lv))$coefficients[2,1:2]
# by = summary(lm(Rt~Lv))$coefficients[2,1:2]
# mr_in = mr_input(bx[1],bx[2],by[1],by[2])
# res = mr_ivw(mr_in)
# res@Pvalue

# Try rcausal
# Installation: 
# install.packages("stringr")
# install.packages("rJava")
# install.packages(devtools)
# library(devtools)
# install_github("bd2kccd/r-causal")
install.packages("rJava",,"http://rforge.net/",type="source")
library(rcausal);library(rJava)
.jinit()
df = discDATA[,inds]
rows = !apply(is.na(df),1,any)
fges_model = fci(df[rows,],continuous = F, verbose = TRUE,significance = 1e-5)
ee = fges_model$edges
ee[grepl("GWS",ee)]

#####################################################################
# Separation or reversal via cond independence examples
# RHR and heart disease via bblockers
summary(lm(RHR~BBLOCKERS+BRADY+TACHY+HeartDISEASE+SEX+AGE,data=DATA))
summary(lm(RHR~BRADY+TACHY+HeartDISEASE+SEX+AGE,data=DATA))
# GWS and heart disease via HEIGHT
summary(glm(AFIB~SEX+AGE+GWS+SMOKING,data=DATA,family = binomial(link='logit')))
summary(glm(HeartDISEASE~SEX+AGE+GWS+HEIGHT+factor(SMOKING),data=DATA,family = binomial(link='logit')))
# GWS scores vs. AFIB
summary(glm(AFIB~GWSHeight+SEX+AGE+HEIGHT,family = binomial(link='logit'),data=DATA))
summary(glm(AFIB~GWSnonheight+SEX+AGE+HEIGHT,family = binomial(link='logit'),data=DATA))
summary(glm(DIABETES~GWSHeight+SEX+AGE,family = binomial(link='logit'),data=DATA))
summary(glm(DIABETES~GWSnonheight+SEX+AGE,family = binomial(link='logit'),data=DATA))
summary(glm(MI_Disease~GWSnonheight+SEX+AGE,family = binomial(link='logit'),data=DATA))
summary(glm(HYPERTENSION~GWSnonheight+GWSHeight+SEX+AGE,family = binomial(link='logit'),data=DATA))
summary(lm(RHR~GWSnonheight+SEX+AGE,data=DATA))
# GWS vs. EHR: positive association
EHR_D = data.frame(EHR,DATA[names(EHR),])
summary(lm(EHR~GWSHeight+GWSnonheight+SEX+AGE,data=EHR_D))
summary(lm(EHR~GWSHeight+GWSnonheight+SEX+AGE+HEIGHT+AFIB,data=EHR_D))
# GWS vs. Height: negative association
summary(lm(HEIGHT~GWSnonheight+GWSHeight+SEX+AGE,data=DATA))
# HEIGHT vs. AFIB: positive association
summary(glm(AFIB~SEX+AGE+HEIGHT,family = binomial(link='logit'),data=DATA))
summary(glm(AFIB~EHR+SEX+AGE,family = binomial(link='logit'),data=EHR_D))
summary(glm(AFIB~EHR+SEX+AGE+HEIGHT,family = binomial(link='logit'),data=EHR_D))
summary(glm(AFIB~EHR+SEX+AGE+HEIGHT+FAT,family = binomial(link='logit'),data=EHR_D))
summary(glm(AFIB~EHR+SEX+AGE+HEIGHT+FAT+DBP+INCOME,family = binomial(link='logit'),data=EHR_D))
summary(lm(HEIGHT~SEX+AGE+AFIB,data=DATA))
##############################################################################

# 
# # Recovery data
# inds = intersect(subjects_with_data,names(expo_vs$Recovery))
# rec_dframe = data.frame(rec=expo_vs$Recovery[inds],age=AGE[inds],sex=SEX[inds],
#                                is_healthy=is_healthy[inds],HeartDISEASE=HeartDISEASE[inds],
#                                fat=FAT[inds],HEIGHT=HEIGHT[inds])
# summary(lm(rec~.,data=rec_dframe))
# summary(lm(rhr~.,data=rhr_dframe[inds,]))
# summary(glm(HeartDISEASE~HEIGHT+SEX))
# boxplot(ehr~BBLOCKERS+sex,data=exercisehr_dframe)
# 
# cor(gws1[inds],expo_vs$ExerciseHR[inds])
# cor(gws2[inds],expo_vs$ExerciseHR[inds])

# Look at correlations of gws with the covariates
cov_corrs_all = apply(covariate_matrix[,names(which(feature_is_numeric))],2,pairwise_cor,y=gws2)
cov_corrs_ehr_subjs = apply(covariate_matrix[names(EHR),names(which(feature_is_numeric))],
                            2,pairwise_cor,y=gws2[names(EHR)])
plot(cov_corrs_all,cov_corrs_ehr_subjs);abline(0,1)
sort(abs(cov_corrs_all))

# # For each disease look at sick vs healthy:
# # Exercise HR, RHR, gws
# get_disease_vs_scores_est<-function(x,y,sampsize=1000){
#   inds = intersect(names(x),names(y))
#   x = x[inds];y=y[inds]
#   x1 = x[y==0]
#   x2 = x[y==1]
#   if(length(x2)<5){return(c(NA,NA,NA))}
#   lmobj = summary(lm(x~y))$coefficients
#   beta = lmobj[2,1]
#   betap = lmobj[2,4]
#   samp1 = sample(x1,min(length(x1),sampsize))
#   samp2=sample(x2,min(length(x2),sampsize))
#   pval = wilcox.test(samp1,samp2)$p.value
#   return(c(beta,betap,pval))
# }
# library(speedglm)
# get_disease_vs_scores_est_with_covs<-function(x,y,covs,useglm=F){
#   inds = intersect(names(x),names(y))
#   inds = intersect(inds,rownames(covs))
#   x = x[inds];y=y[inds];covs=covs[inds,]
#   d=data.frame(x=x,y=y,covs=covs)
#   if(!useglm){
#     lmobj = summary(lm(x~.,data=d))$coefficients
#     betay = lmobj["y",1]
#     betayp = lmobj["y",4]
#     return(c(betay,betayp))
#   }
#   logiobj = speedglm(factor(y)~.,family=binomial(link='logit'),data=d)
#   logiobj = summary(logiobj)
#   logiobj = as.matrix(logiobj$coefficients)
#   mode(logiobj) = 'numeric'
#   betax = logiobj["x",1]
#   betaxp = logiobj["x",4]
#   return(c(betax,betaxp))
# }
# icd_cols_corr_analysis = c()
# for(j in 1:ncol(icd_data)){
#   disease_v = icd_data[,j]; names(disease_v) = rownames(icd_data)
#   expo_res = get_disease_vs_scores_est(expo_vs$ExerciseHR,disease_v)
#   gws_res = get_disease_vs_scores_est(gws,disease_v)
#   rhr_res = get_disease_vs_scores_est(rhr,disease_v)
#   icd_cols_corr_analysis = rbind(icd_cols_corr_analysis,
#                                  c(expo_res,gws_res,rhr_res))
#   rownames(icd_cols_corr_analysis)[nrow(icd_cols_corr_analysis)] = colnames(icd_data)[j]
#   print(colnames(icd_data)[j])
#   print(rhr_res)
# }
# ord = order(icd_cols_corr_analysis[,9])
# icd_cols_corr_analysis[ord,]
# ord = order(icd_cols_corr_analysis[,6])
# icd_cols_corr_analysis[ord,]
# 
# icd_cols_adjusted_corr_analysis = c()
# batch_inds = grepl(colnames(simple_covs),pattern="batch")
# for(j in 1:ncol(icd_data)){
#   disease_v = icd_data[,j]
#   names(disease_v) = rownames(icd_data)
#   expo_res = get_disease_vs_scores_est_with_covs(expo_vs$ExerciseHR,disease_v,simple_covs[,!batch_inds])
#   gws_res = get_disease_vs_scores_est_with_covs(gws,disease_v,simple_covs[,!batch_inds])
#   rhr_res = get_disease_vs_scores_est_with_covs(rhr,disease_v,simple_covs[,!batch_inds])
#   icd_cols_adjusted_corr_analysis = rbind(icd_cols_adjusted_corr_analysis,
#                                           c(expo_res,gws_res,rhr_res))
#   rownames(icd_cols_adjusted_corr_analysis)[nrow(icd_cols_adjusted_corr_analysis)] = colnames(icd_data)[j]
#   print(colnames(icd_data)[j])
#   print(rhr_res)
# }
# ord = order(icd_cols_adjusted_corr_analysis[,6])
# icd_cols_adjusted_corr_analysis[ord,]
# ord = order(icd_cols_corr_analysis[,6])
# icd_cols_corr_analysis[ord,]
# selected_inds = icd_cols_adjusted_corr_analysis[,4]<1e-5
# table(selected_inds)
# icd_cols_adjusted_corr_analysis[selected_inds,]
# table(icd_cols_adjusted_corr_analysis[,5]>0, icd_cols_adjusted_corr_analysis[,6]<0.01)

