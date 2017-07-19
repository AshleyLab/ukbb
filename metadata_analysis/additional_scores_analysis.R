## NOTE!!!
# This script is under development and has many parts copied from
# phenotype_and_covariate_scores_analysis.R

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

# define constants
MAX_NUM_CATEGORIES_IN_DISCRETE_COVARIATE = 100
MIN_CLASS_SIZE_IN_DISCRETE_COVARIATE = 500

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
###############################################
################# Phenotypes ##################
###############################################
###############################################
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
status_cols = get_regex_cols(colnames(pheno_data),"status of test")
additional_scores[['exercise_test_status']] = pheno_data[,status_cols]
save(additional_scores,file = "metadata_analysis_additional_scores")
gc()

# For each additional score: simplify, dummy vars, impute, and run PCA
additional_scores_pcs = list()
for(nn in names(additional_scores)){
  mat = additional_scores[[nn]]
  all_nas = apply(is.na(mat),1,all)
  mat = mat[!all_nas,]
  print(dim(mat))
  print(colnames(mat))
  mat = merge_columns_by_their_features(mat,feature_code2name[colnames(mat)])
  if(is.null(dim(mat))){
    additional_scores_pcs[[nn]] = mat
    next
  }
  print(dim(mat))
  print(colnames(mat))
  col_percent_na = colSums(is.na(mat))/nrow(mat)
  mat = mat[,col_percent_na<0.5]
  print(dim(mat))
  print(colnames(mat))
  all_nas = apply(is.na(mat),1,all)
  mat = mat[!all_nas,]
  print(dim(mat))
  print(colnames(mat))
  initial_pca_success = F
  try({
    pcs = get_pcs_for_score_set(mat,num_pcs=2)
    initial_pca_success = T
  })
  if(!initial_pca_success){
    has_na = rowSums(is.na(mat))>0
    mat = mat[!has_na,]
    pcs = get_pcs_for_score_set(mat,num_pcs=2)
  }
  colnames(pcs) = paste(nn,"_PC",1:ncol(pcs),sep="")
  additional_scores_pcs[[nn]] = list(pcs=pcs,initial_pca_success=initial_pca_success)
}
save(additional_scores,additional_scores_pcs,file = "metadata_analysis_additional_scores")

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
subject_set = unique(unlist(sapply(additional_scores,rownames)))

# Filter 1: Define data column types to be excluded automatically
categories_to_exclude = c("Cognitive function")
potential_subcategories_to_exclude = c()
subcategories_to_exclude = c(
  "ECG during exercise","Raw accelerometer statistics","Accelerometer wear time duration","Genotype calls &amp; imputation",
  "Brain MRI","Acceleration averages"
)
# Add the columns of the additional PA scores above
# Filter 1
features_to_exclude = unique(unlist(sapply(additional_scores,colnames)))
# Filter 2
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
print(dim(pheno_data))
gc()

# Filter 1.1: Features with 100% NAs in the resulting matrix (save some running time)
percent_nas_col = colSums(is.na(pheno_data))/nrow(pheno_data)
pheno_data = pheno_data[,percent_nas_col < 1]
print(dim(pheno_data))
gc()

pheno_data_feature2name = feature_code2name[colnames(pheno_data)]
colnames(pheno_data)[is.na(pheno_data_feature2name)]
# Corrections for the NAs in the names
names(pheno_data_feature2name) = colnames(pheno_data)
pheno_data_feature2name["Weight.0.0.1"] = "Weight"
pheno_data_feature2name["Weight.1.0.1"] = "Weight"
pheno_data_feature2name["Body mass index (BMI).0.0.1"] = "Body mass index (BMI)"
pheno_data_feature2name["Body mass index (BMI).1.0.1"] = "Body mass index (BMI)"

# Filter 6: exclude categorical features with too many classes (>100) 
# or that all classes are too small (<500)
feature_is_numeric = c()
for(j in 1:ncol(pheno_data)){
  feature_is_numeric[j] = !all(is.na(as.numeric(as.character(pheno_data[,j]))))
}
print(table(feature_is_numeric))
names(feature_is_numeric) = colnames(pheno_data)
features_to_exclude = c()
for(j in 1:ncol(pheno_data)){
  if(feature_is_numeric[j]){next}
  if(length(unique(pheno_data[,j]))>MAX_NUM_CATEGORIES_IN_DISCRETE_COVARIATE){
    features_to_exclude = c(features_to_exclude,j)
    next
  }
  v = pheno_data[,j];v = v[!is.na(v)]
  if(max(table(v))<MIN_CLASS_SIZE_IN_DISCRETE_COVARIATE){
    features_to_exclude = c(features_to_exclude,j)
  }
}
length(features_to_exclude)
colnames(pheno_data)[features_to_exclude]
pheno_data = pheno_data[,-features_to_exclude]
pheno_data_feature2name = pheno_data_feature2name[colnames(pheno_data)]
print(dim(pheno_data))
gc()

# Filter 3+4: Merge by feature name, exclude certain regular expressions
# The merge is done by taking the first column and using the others to fill NAs
regex_to_exclude = c("device","last hour","interpolated","operation code",
                     "method of recording","genetic principal","quality control","date of")
new_pheno_dat = data.frame(subjnames = rownames(pheno_data))
for(nn in unique(pheno_data_feature2name)){
  if(any(sapply(regex_to_exclude,grepl,x=nn,ignore.case=T))){next}
  curr_cols = pheno_data[,pheno_data_feature2name == nn]
  if(sum(pheno_data_feature2name == nn)==1){
    v = data.frame(nn = curr_cols)
  }
  else{
    v = data.frame(nn = merge_feature_columns(curr_cols))
  }
  new_pheno_dat = cbind(new_pheno_dat,v)
  colnames(new_pheno_dat)[ncol(new_pheno_dat)] = nn
  #print(ncol(new_pheno_dat))
}
rownames(new_pheno_dat) = rownames(pheno_data)
new_pheno_dat = new_pheno_dat[,-1]
print(dim(new_pheno_dat))
gc()

# Filter 5: For covariate analysis, exclude features that are 20% NAs or higher
features_to_exclude = apply(is.na(new_pheno_dat),2,sum)/nrow(new_pheno_dat) >= 0.2
new_pheno_dat = new_pheno_dat[,!features_to_exclude]
dim(new_pheno_dat)
gc()

covariate_matrix = new_pheno_dat
feature_is_numeric = c(); num_categories = c()
for(j in 1:ncol(covariate_matrix)){
  feature_is_numeric[j] = !all(is.na(as.numeric(as.character(covariate_matrix[,j]))))
  if(feature_is_numeric[j]){num_categories[j] = NA}
  else{num_categories[j] = length(unique(covariate_matrix[,j]))}
}
table(feature_is_numeric)
table(num_categories)
names(feature_is_numeric) = colnames(covariate_matrix)
save(covariate_matrix,feature_is_numeric,file="covariate_matrix.RData")
rm(pheno_data,new_pheno_dat)
gc()

# Check that sex, age, height, weight are there
get_regex_cols(colnames(covariate_matrix),"age",ignore.case=T)
get_regex_cols(colnames(covariate_matrix),"sex",ignore.case=T)
get_regex_cols(colnames(covariate_matrix),"height",ignore.case=T)
get_regex_cols(colnames(covariate_matrix),"weight",ignore.case=T)
get_regex_cols(colnames(covariate_matrix),"bmi",ignore.case=T)
tmp_feature = get_regex_cols(colnames(covariate_matrix),"smok",ignore.case=T)[2]

# # 1.1 SText 2 (The data)
# # map to pheno_data and look at the raw features vs the computed ones
# tmp_feature = get_regex_cols(colnames(covariate_matrix),"smok",ignore.case=T)[2]
# tmp_raw_features = names(which(feature_code2name == tmp_feature))
# feature_subcategory_data[tmp_feature]
# tmp_raw_is_num = !all(is.na(pheno_data[,tmp_raw_features[1]]))
# tmp_raw_non_na_sets = apply(!is.na(pheno_data[,tmp_raw_features]),2,which)
# sapply(tmp_raw_non_na_sets,length)
# covered_samples = unique(unlist(tmp_raw_non_na_sets))
# length(covered_samples)
# raw_feature_tables = apply(pheno_data[,tmp_raw_features],2,table)
# rowSums(raw_feature_tables)
# # Look at different merging options
# merged_col1 = merge_feature_columns(pheno_data[,tmp_raw_features])
# table(merged_col1)
# merged_col2 = merge_feature_columns(pheno_data[,tmp_raw_features],take_first=F)
# table(merged_col2)
# merged_col3 = apply(pheno_data[,tmp_raw_features],1,max,na.rm=T)
# table(merged_col3)
# merged_col4 = apply(pheno_data[,tmp_raw_features],1,min,na.rm=T)
# table(merged_col4)
# merged_col5 = apply(pheno_data[,tmp_raw_features],1,mean,na.rm=T)
# table(merged_col5)
# ###
# # check consistency:
# inter1 = intersect(tmp_raw_na_sets[[1]],tmp_raw_na_sets[[2]])
# inter1_matrix = pheno_data[inter1,tmp_raw_features[1:2]]
# table(apply(inter1_matrix,1,function(x)(all(x==x[1]))))
# length(intersect(tmp_raw_na_sets[[1]],tmp_raw_na_sets[[3]]))
# computed_tmp_feature = covariate_matrix[,tmp_feature[5]]


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

# Conservative analysis
fitness_vs_covs_lm_objects = list()
for (fitness_score_ind in 1:ncol(fitness_scores_matrix)){
  y = fitness_scores_matrix[,fitness_score_ind]
  currname = colnames(fitness_scores_matrix)[fitness_score_ind]
  lm1 = get_lm_residuals(y,covariate_matrix,use_categorical=T,max_num_classes=20,feature_is_numeric=feature_is_numeric)
  plot(y[names(lm1[[1]]$residuals)],lm1[[1]]$residuals,ylab="Residual",xlab="Fitness score")
  fitness_vs_covs_lm_objects[[currname]] = lm1
  save(fitness_vs_covs_lm_objects,file="fitness_analysis_fitness_vs_covs_lm_objects.RData")
}
r2_scores = sapply(fitness_vs_covs_lm_objects,function(x)summary(x$lm)[["r.squared"]])

# Plot for the report: Figure S2.2 
par(mfrow=c(2,2))
r2_scores = sapply(fitness_vs_covs_lm_objects,function(x)summary(x$lm)[["r.squared"]])
r2_scores = format(r2_scores,digits=2)
for (fitness_score_ind in 1:ncol(fitness_scores_matrix)){
  y = fitness_scores_matrix[,fitness_score_ind]
  currname = colnames(fitness_scores_matrix)[fitness_score_ind]
  lm1 = fitness_vs_covs_lm_objects[[currname]]
  plot(y[names(lm1[[1]]$residuals)],lm1[[1]]$residuals,ylab="Residual",xlab="Fitness score",
       main=paste(currname,", R^2=",r2_scores[fitness_score_ind],sep=""))
}

accelerometry_scores_to_residuals = list()
for (j in 1:ncol(accelerometry_scores)){
  y = accelerometry_scores[,j]
  names(y) = rownames(accelerometry_scores)
  currname = colnames(accelerometry_scores)[j]
  lm1 = get_lm_residuals(y,covariate_matrix,use_categorical=T,max_num_classes=20,feature_is_numeric=feature_is_numeric)
  plot(y[names(lm1[[1]]$residuals)],lm1[[1]]$residuals,ylab="Residual",xlab="Fitness score")
  accelerometry_scores_to_residuals[[currname]] = lm1[[1]]$residuals
  save(accelerometry_scores_to_residuals,file="accelereometry_analysis_score_vs_covs_residuals_conservative.RData")
}
# TBD
additional_scores_residuals = list(){}

# Simple analysis
simple_covs = covariate_matrix[,c("Sex","Age when attended assessment centre","Body mass index (BMI)","BMI","Standing height")]
fitness_vs_covs_lm_objects = list()
for (fitness_score_ind in 1:ncol(fitness_scores_matrix)){
  y = fitness_scores_matrix[,fitness_score_ind]
  currname = colnames(fitness_scores_matrix)[fitness_score_ind]
  lm1 = get_lm_residuals(y,simple_covs,use_categorical=T,max_num_classes=20,feature_is_numeric=feature_is_numeric)
  plot(y[names(lm1[[1]]$residuals)],lm1[[1]]$residuals,ylab="Residual",xlab="Fitness score")
  fitness_vs_covs_lm_objects[[currname]] = lm1
  save(fitness_vs_covs_lm_objects,file="fitness_analysis_fitness_vs_covs_lm_objects_simple.RData")
}
accelerometry_scores_to_residuals = list()
for (j in 1:ncol(accelerometry_scores)){
  y = accelerometry_scores[,j]
  names(y) = rownames(accelerometry_scores)
  currname = colnames(accelerometry_scores)[j]
  lm1 = get_lm_residuals(y,simple_covs,use_categorical=T,max_num_classes=20,feature_is_numeric=feature_is_numeric)
  plot(y[names(lm1[[1]]$residuals)],lm1[[1]]$residuals,ylab="Residual",xlab="Fitness score")
  accelerometry_scores_to_residuals[[currname]] = lm1[[1]]$residuals
  save(accelerometry_scores_to_residuals,file="accelereometry_analysis_score_vs_covs_residuals_simple.RData")
}
# TBD
additional_scores_residuals = list(){}

# Plot for the report: Figure S2.2 
load("fitness_analysis_fitness_vs_covs_lm_objects_simple.RData")
par(mfrow=c(2,2))
r2_scores = sapply(fitness_vs_covs_lm_objects,function(x)summary(x$lm)[["r.squared"]])
r2_scores = format(r2_scores,digits=2)
for (fitness_score_ind in 1:ncol(fitness_scores_matrix)){
  y = fitness_scores_matrix[,fitness_score_ind]
  currname = colnames(fitness_scores_matrix)[fitness_score_ind]
  lm1 = fitness_vs_covs_lm_objects[[currname]]
  plot(y[names(lm1[[1]]$residuals)],lm1[[1]]$residuals,ylab="Residual",xlab="Fitness score",
       main=paste(currname,", R^2=",r2_scores[fitness_score_ind],sep=""))
}




