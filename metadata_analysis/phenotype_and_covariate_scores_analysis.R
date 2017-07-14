try({setwd("/Users/davidhsu/Documents/ukbb")})
try({setwd("/Users/david/Desktop/ukbb/")})
try({setwd('/scratch/PI/euan/projects/ukbb/da_dh/')})
# We cannot change the working directory as we work in RStudio, so we comment out this line temporarily.
source("auxiliary_functions.R")
library(xlsx)

###############################################
###############################################
############## Load feature data ##############
###############################################
###############################################

# Load the preprocessed data and the pheno data
load("biobank_collated_pheno_data.RData")

# Analysis of the features
feature_code2name = sapply(colnames(pheno_data),function(x)strsplit(x,split='\\.')[[1]][1])

# Load mapping of pheno data column names into categories
feature_metadata =  read.csv('Data_Dictionary_Showcase.csv')
feature_subcategory_data = as.character(feature_metadata[,"Cat3_Title"])
names(feature_subcategory_data) = as.character(feature_metadata[,"Field"])
feature_category_data = as.character(feature_metadata[,"Cat2_Title"])
names(feature_category_data) = as.character(feature_metadata[,"Field"])

# Load the manual analysis of the feature names, done by Daryl
feature_data_recommentadions = read.xlsx2("HR_exercise_pred_100_covariate_analysis.xlsx",1)
feature_recommentadions = as.character(feature_data_recommentadions[,ncol(feature_data_recommentadions)])
names(feature_recommentadions) = feature_data_recommentadions[,1]
feature_recommentadions = feature_recommentadions[names(feature_recommentadions)!=""]
features_to_exclude_manual_analysis = names(feature_recommentadions)[grepl("ignore",feature_recommentadions,ignore.case = T)
                                                                     | grepl("score",feature_recommentadions,ignore.case = T)]
features_to_exclude_manual_analysis = unname(sapply(features_to_exclude_manual_analysis,function(x)strsplit(x,split='\\.')[[1]][1]))

###############################################
#################### End ######################
###############################################
###############################################

###############################################
###############################################
################# Phenotypes ##################
###############################################
###############################################

# Fitness scores:
load("fitness_analysis_final_fitness_scores.RData")
# Additional phenotypes:
# a. Anna's scores
accelerometry_scores = read.delim('accelerometry_phenotypes_anna/accelerometery_continuous_features_quantile_normalized_no_outliers.tsv')
all(accelerometry_scores[,1]==accelerometry_scores[,2])
rownames(accelerometry_scores) = accelerometry_scores[,1]
accelerometry_scores = accelerometry_scores[,-c(1:2)]
accelerometry_scores_all_nas = apply(is.na(accelerometry_scores),1,all)
table(accelerometry_scores_all_nas)
accelerometry_scores = accelerometry_scores[!accelerometry_scores_all_nas,]
# b. IDC-9 (???)
# c. PA related columns
extract_feature_cols_by_category<-function(cols,cols2names,names2cats,category){
  curr_features = c()
  for(nn in cols){
    curr_name = cols2names[nn]
    curr_cat = names2cats[curr_name]
    if(is.na(curr_cat)){next}
    if(curr_cat==category){
      curr_features = c(curr_features,nn)
    }
  }
  return(curr_features)
}
additional_scores = list()
spiro_cols = extract_feature_cols_by_category(colnames(pheno_data),feature_code2name,feature_subcategory_data,"Spirometry")
additional_scores[['Spirometry']] = pheno_data[,spiro_cols]
dxa_cols = extract_feature_cols_by_category(colnames(pheno_data),feature_code2name,feature_subcategory_data,"DXA assessment")
additional_scores[['DXA']] = pheno_data[,dxa_cols]
impedence_cols = get_regex_cols(colnames(pheno_data),"Impedance",ignore.case=T)
feature_category_data[feature_code2name[impedence_cols]]
additional_scores[['Impedance']] = pheno_data[,impedence_cols]
pulse_rate_cols = get_regex_cols(colnames(pheno_data),"pulse rate\\.",ignore.case=T)
additional_scores[['Pulse rate']] = pheno_data[,pulse_rate_cols]
hand_grip_cols = get_regex_cols(colnames(pheno_data),"Hand grip",ignore.case=T)
additional_scores[['hand grip']] = pheno_data[,hand_grip_cols]
gc()

###############################################
#################### End ######################
###############################################
###############################################

###############################################
###############################################
############## Filter features ################
###############################################
###############################################

# Define the sample set to work with: union of the above
intersect(rownames(accelerometry_scores),rownames(fitness_scores_matrix))
subject_set = union(rownames(accelerometry_scores),rownames(fitness_scores_matrix))

# Define data column types to be excluded automatically
categories_to_exclude = c("Cognitive function")
potential_subcategories_to_exclude = c()
subcategories_to_exclude = c(
  "ECG during exercise","Raw accelerometer statistics","Accelerometer wear time duration","Genotype calls &amp; imputation",
  "Brain MRI","Acceleration averages"
)
# Add the columns of the additional PA scores above
features_to_exclude = unique(unlist(sapply(additional_scores,colnames)))

# Remove columns using the above recommentations and the defined categories to remove
for(j in 1:ncol(pheno_data)){
  currf = colnames(pheno_data)[j]
  currf_name = feature_code2name[currf]
  curr_cat = feature_category_data[currf_name]
  curr_subcat = feature_subcategory_data[currf_name]
  if(is.element(curr_cat,set=categories_to_exclude)){
    features_to_exclude = c(features_to_exclude,currf)
    next
  }
  if(is.element(curr_subcat,set=subcategories_to_exclude)){
    features_to_exclude = c(features_to_exclude,currf)
    next
  }
  if(is.element(currf_name,set=features_to_exclude_manual_analysis)){
    features_to_exclude = c(features_to_exclude,currf)
    next
  }
}
features_to_exclude = unique(features_to_exclude)
length(features_to_exclude)

# Filter 1: Trim the data matrix by columns and rows above
pheno_data = pheno_data[,-which(is.element(colnames(pheno_data),set=features_to_exclude))]
pheno_data = pheno_data[subject_set,]
gc()
# Filter 2: Features with >99% NAs in the resulting matrix
percent_nas_col = colSums(is.na(pheno_data))/nrow(pheno_data)
pheno_data = pheno_data[,percent_nas_col < 0.99]
pheno_data_feature2name = feature_code2name[colnames(pheno_data)]
# Corrections for the NAs in the names
names(pheno_data_feature2name) = colnames(pheno_data)
pheno_data_feature2name["Weight.0.0.1"] = "Weight"
pheno_data_feature2name["Weight.1.0.1"] = "Weight"
pheno_data_feature2name["Body mass index (BMI).0.0.1"] = "BMI"
pheno_data_feature2name["Body mass index (BMI).1.0.1"] = "BMI"
gc()

# Filter 3: exclude categorical features with too many classes (>100) 
# or that all classes are too small (<500)
feature_is_numeric = c()
for(j in 1:ncol(pheno_data)){
  feature_is_numeric[j] = !all(is.na(as.numeric(as.character(pheno_data[,j]))))
}
table(feature_is_numeric)
names(feature_is_numeric) = colnames(pheno_data)
features_to_exclude = c()
for(j in 1:ncol(pheno_data)){
  if(feature_is_numeric[j]){next}
  if(length(unique(pheno_data[,j]))>100){
    features_to_exclude = c(features_to_exclude,j)
    next
  }
  v = pheno_data[,j]
  v = v[!is.na(v)]
  if(max(table(v))<500){
    features_to_exclude = c(features_to_exclude,j)
  }
}
length(features_to_exclude)
pheno_data = pheno_data[,-features_to_exclude]
pheno_data_feature2name = pheno_data_feature2name[colnames(pheno_data)]

# Filter 4: Merge by feature name, exclude certain regular expressions
# The merge is done by taking the first column and using the others to fill NAs
regex_to_exclude = c("device","last hour","interpolated","operation code",
                     "method of recording","genetic principal","quality control","date of")
new_pheno_dat = c()
for(nn in unique(pheno_data_feature2name)){
  if(any(sapply(regex_to_exclude,grepl,x=nn,ignore.case=T))){next}
  curr_cols = pheno_data[,pheno_data_feature2name == nn]
  if(sum(pheno_data_feature2name == nn)==1){
    new_pheno_dat = cbind(new_pheno_dat,curr_cols)
  }
  else{
    v = curr_cols[,1]
    for(j in 2:ncol(curr_cols)){
      curr_NAs = is.na(v)
      v[curr_NAs] = curr_cols[curr_NAs,j]
    }
    new_pheno_dat = cbind(new_pheno_dat,v)
  }
  colnames(new_pheno_dat)[ncol(new_pheno_dat)] = nn
}
dim(new_pheno_dat)
rownames(new_pheno_dat) = rownames(pheno_data)

# Filter 5: For covariate analysis, exclude features that are 40%
features_to_exclude = apply(is.na(new_pheno_dat),2,sum)/nrow(new_pheno_dat) >= 0.2
new_pheno_dat = new_pheno_dat[,!features_to_exclude]
dim(new_pheno_dat)
# Check that sex, age, height, weight are there
get_regex_cols(colnames(new_pheno_dat),"age",ignore.case=T)
get_regex_cols(colnames(new_pheno_dat),"sex",ignore.case=T)
get_regex_cols(colnames(new_pheno_dat),"height",ignore.case=T)
get_regex_cols(colnames(new_pheno_dat),"weight",ignore.case=T)
get_regex_cols(colnames(new_pheno_dat),"bmi",ignore.case=T)
get_regex_cols(colnames(new_pheno_dat),"smok",ignore.case=T)

covariate_matrix = new_pheno_dat
feature_is_numeric = c()
for(j in 1:ncol(covariate_matrix)){
  feature_is_numeric[j] = !all(is.na(as.numeric(as.character(covariate_matrix[,j]))))
}
table(feature_is_numeric)
names(feature_is_numeric) = colnames(covariate_matrix)

save(covariate_matrix,feature_is_numeric,file="covariate_matrix.RData")

rm(pheno_data,new_pheno_dat)
gc()

###############################################
###############################################
#################### End ######################
###############################################
###############################################

###############################################
###############################################
################## Residuals ##################
###############################################
###############################################

fitness_vs_covs_lm_objects = list()
for (fitness_score_ind in 1:ncol(fitness_scores_matrix)){
  y = fitness_scores_matrix[,fitness_score_ind]
  currname = colnames(fitness_scores_matrix)[fitness_score_ind]
  if(is.element(currname,set=names(fitness_vs_covs_lm_objects))){next}
  # All covariates, without the additionally excluded
  lm1 = get_lm_residuals(y,covariate_matrix,use_categorical=T,max_num_classes=20,feature_is_numeric=feature_is_numeric)
  plot(y[names(lm1[[1]]$residuals)],lm1[[1]]$residuals,ylab="Residual",xlab="Fitness score")
  lm1_summ = summary(lm1[[1]])
  sigs = lm1_summ$coefficients[,4]
  sort(sigs)[1:5]
  fitness_vs_covs_lm_objects[[currname]] = lm1
  save(fitness_vs_covs_lm_objects,file="fitness_analysis_fitness_vs_covs_lm_objects.RData")
}

# Load the genetic PCA data and print the tables for the GWAS
# TODO:
# 1. Discretize?
# 2. What to do with NAs?
genetic_pca_data_path = "plink/may16.eigenvec"
genetic_pca_data = read.delim(genetic_pca_data_path,sep=" ",header = F)
dim(genetic_pca_data)
gwas_data_samples = as.character(genetic_pca_data[,2])
gwas_data_pcs = genetic_pca_data[,3:4]
colnames(gwas_data_pcs) = paste("PC",1:ncol(gwas_data_pcs),sep="")
gwas_data_residuals = c()
# this vector specifies the direction of the scores - whether higher is 
# better fitness or not
for(nn in names(fitness_vs_covs_lm_objects)){
  lm_obj = fitness_vs_covs_lm_objects[[nn]][[1]]
  curr_res = lm_obj$residuals
  curr_res = curr_res[gwas_data_samples]
  names(curr_res) = gwas_data_samples
  NA_samples = gwas_data_samples[is.na(curr_res)]
  gwas_data_residuals = cbind(gwas_data_residuals,curr_res)
}
colnames(gwas_data_residuals) = paste("Residuals_",names(fitness_vs_covs_lm_objects),sep="")
colnames(gwas_data_residuals) = gsub(colnames(gwas_data_residuals),pattern=" ",replace="_")
gwas_data = cbind(gwas_data_residuals,as.matrix(gwas_data_pcs))
get_pairwise_corrs(gwas_data)
save(gwas_data,file="June14_2017_gwas_data_table.RData")
write.table(gwas_data,file = "June14_2017_fitness_scores_gwas_data_table.txt",sep="\t",quote=F)

###############################################
###############################################
#################### End ######################
###############################################
###############################################

###############################################
###############################################
########## Fitness scores analysis ############
###############################################
###############################################

# TODO: check later
library(corrplot)
par(mfrow=c(2,2),mar=c(2,4,2,2))
boxplot(get_list_by_values(subject_ols_preds_100,subject_technical_class,1,"No data"),main="HR_pred 100")
boxplot(get_list_by_values(rest_HR_decline_score[,2],subject_technical_class,1,"No data"),main="Rest_diff")
boxplot(get_list_by_values(rest_HR_decline_WD_normalized_score[,2],subject_technical_class,1,"No data"),main="Rest_diff_WD")
boxplot(get_list_by_values(rest_HR_decline_HR_normalized_score[,2],subject_technical_class,1,"No data"),main="Rest_ratio")

# Fitness scores vs. other covariates - check for associations
cut_by_quantiles<-function(x,nbreaks=5){
  if(length(unique(x))<=nbreaks){return(x)}
  probs = seq(0,1,length.out = nbreaks+1)
  qs = unique(quantile(x,probs = probs,na.rm=T))
  xv = cut(x,breaks=qs)
  return(xv)
}

cut_size = 5
Mc = covariate_matrix
for(j in 1:ncol(Mc)){
  if(feature_is_numeric[colnames(Mc)[j]]){
    v = as.numeric(as.character(Mc[,j]))
    v = cut_by_quantiles(v,nbreaks=min(cut_size,length(unique(v))))
    Mc[,j] = v
  }
}

library(entropy)
covariance_correlation_summary_tables = list()
curr_subject_set = intersect(rownames(covariate_matrix),rownames(fitness_scores_matrix))
for(fit_ind in 1:ncol(fitness_scores_matrix)){
  fv = fitness_scores_matrix[curr_subject_set,fit_ind]
  fv_disc = cut_by_quantiles(fv,nbreaks = cut_size)
  fv_is_na = is.na(fv)
  summary_table = c()
  for(j in 1:ncol(Mc)){
    covv = covariate_matrix[curr_subject_set,j]
    covv_disc = Mc[curr_subject_set,j]
    covv_na = is.na(covv)
    cov_name = colnames(Mc)[j]
    cov_category = feature_category_data[cov_name]
    # correlations
    non_na_inds = !fv_is_na & !is.na(covv)
    sp_rho = NA ; sp_rho_p = NA
    if(sum(non_na_inds)>1000 & feature_is_numeric[cov_name]){
      x1 = as.numeric(covv[non_na_inds])
      x2 = fv[non_na_inds]
      sp_rho = cor(x1,x2,method='spearman')
      sp_rho_p = cor.test(x1,x2,method='spearman')$p.value
    }
    # na vs na
    na_na_table = table(fv_is_na,covv_na)
    na_vs_na_chisq_pvalue = NA;na_vs_na_MI = NA
    if(length(na_na_table)==4){
      na_vs_na_chisq_pvalue = chisq.test(na_na_table)$p.value
      na_vs_na_MI = mi.empirical(na_na_table)
    }
    # discrete vs discrete
    disc_chisq_pvalue = NA;disc_MI = NA
    if(length(unique(covv_disc))<100){
      tab = table(as.character(covv_disc),fv_disc)
      disc_chisq_pvalue = chisq.test(tab)$p.value
      disc_MI = mi.empirical(tab)
    }
    # fv na vs. non na scores
    non_na_cov_inds = !covv_na
    fv_na_chisq_p = NA; fv_na_mi=NA
    if(sum(non_na_cov_inds)>100 && length(table(fv_is_na))>1 && length(unique(covv_disc))<100){
      na_tab = table(as.character(covv_disc[non_na_cov_inds]),fv_is_na[non_na_cov_inds])
      fv_na_chisq_p = chisq.test(na_tab)$p.value
      fv_na_mi = mi.empirical(na_tab)
    }
    
    cov_summary_scores = c(disc_MI,disc_chisq_pvalue,
                           fv_na_mi,fv_na_chisq_p,
                           sp_rho,sp_rho_p,
                           na_vs_na_MI,na_vs_na_chisq_pvalue)
    summary_table = rbind(summary_table,c(cov_name,cov_category,cov_summary_scores))
  }
  covariance_correlation_summary_tables[[colnames(fitness_scores)[fit_ind]]] = summary_table
  break
}
names(covariance_correlation_summary_tables) = colnames(fitness_scores_matrix)
for(nn in corrected_names){
  colnames(covariance_correlation_summary_tables[[nn]]) = c("Feature","Category","MI-discretized","ChisqP-discretized",
      "MI-NA fitness","ChisqP-NA fitness","Spearman rho","Spearman rho p","MI-NAs", "ChisqP-NAs")
  write.table(covariance_correlation_summary_tables[[nn]],file=paste(nn,"_cov_analysis.txt",sep=""),sep="\t",quote=F,row.names=F)
}
save(covariance_correlation_summary_tables,file="fitness_analysis_covariance_correlation_summary_tables.RData")














par(mfrow=c(1,1),mar=c(5,5,5,5))
plot(x=pmin(500,-log(na_pval_mat[2,])),y=pmin(500,-log(chisq_pval_mat[2,])),
     xlab="HR_pred_100 - NA (-log 10 p-values)",ylab="HR_pred_100 - scores (-log 10 p-values)",pch=3,lwd=1.5);abline(0,1)

categories = initial_cov_matrix_feature_subcategory
categories[grepl(names(categories),pattern="device",ignore.case=T)] = "device"
pval_list = get_list_by_values(chisq_pval_mat[2,],categories)
pval_list = lapply(pval_list,function(x)x=-log(x,base=10)) 
pval_list = lapply(pval_list,function(x)x[!is.na(x)&!is.nan(x)])
treat_zero_pvals<-function(x){
  x[is.infinite(x)]=300
  return(x)
}
pval_list = lapply(pval_list,treat_zero_pvals)
par(mfrow=c(1,2),mar=c(3,12,3,3))
boxplot(pval_list,horizontal = T,las=2,col="blue",main="HR_pred_100 - scores (-log10 p)")

pval_list = get_list_by_values(na_pval_mat[2,],categories)
pval_list = lapply(pval_list,function(x)x=-log(x,base=10)) 
pval_list = lapply(pval_list,function(x)x[!is.na(x)&!is.nan(x)])
treat_zero_pvals<-function(x){
  x[is.infinite(x)]=300
  return(x)
}
pval_list = lapply(pval_list,treat_zero_pvals)
boxplot(pval_list,horizontal = T,las=2,col="red",main="HR_pred_100 - NA (-log10 p)")

# Some funny plots
tmp = fitness_scores[,2]
ll = get_list_by_values(tmp,drive_faster)
boxplot(ll[c("Never/rarely","Sometimes","Often","Most of the time")],las=2,ylab="Poor fitness",ylim=c(50,250),col=gray.colors(4),main="Fitness vs. driving fast (p < 1E-300)")

sex = as.character(initial_cov_matrix[,"Sex.0.0"])
sex_and_driving = paste(drive_faster,sex,sep=';')
ll2 = get_list_by_values(tmp,sex_and_driving)
strat_classes = c("Never/rarely","Sometimes","Often","Most of the time")
strat_males = paste(strat_classes,"Male",sep=';')
strat_females = paste(strat_classes,"Female",sep=';')
par(mfrow=c(1,2))
boxplot(ll2[strat_males],las=2,ylab="Poor fitness",ylim=c(60,200),col=gray.colors(4),names = strat_classes,main="Males")
boxplot(ll2[strat_females],las=2,ylab="Poor fitness",ylim=c(60,200),col=gray.colors(4),names = strat_classes,main="Females")
