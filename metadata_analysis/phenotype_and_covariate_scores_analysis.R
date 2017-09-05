try({setwd("/Users/davidhsu/Documents/ukbb")})
try({setwd("/Users/david/Desktop/ukbb/")})

# We cannot change the working directory as we work in RStudio, so we comment out this line temporarily.
source("auxiliary_functions.R")
library(xlsx)

gaus_norm<-function(x){
  x_r = (rank(x)-0.5)/length(x)
  x_n = qnorm(x_r)
  return(x_n)
}

###############################################
###############################################
############## Load feature data ##############
###############################################
###############################################

# define constants
MAX_NUM_CATEGORIES_IN_DISCRETE_COVARIATE = 50
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
# b. define a set of additional scores
additional_scores = list()
spiro_cols = extract_feature_cols_by_category(colnames(pheno_data),feature_code2name,feature_subcategory_data,"Spirometry")
additional_scores[['Spirometry']] = pheno_data[,spiro_cols]
dxa_cols = extract_feature_cols_by_category(colnames(pheno_data),feature_code2name,feature_subcategory_data,"DXA assessment")
additional_scores[['DXA']] = pheno_data[,dxa_cols]
impedence_cols = get_regex_cols(colnames(pheno_data),"Impedance",ignore.case=T)
feature_category_data[feature_code2name[impedence_cols]]
additional_scores[['Impedance']] = pheno_data[,impedence_cols]
pulse_rate_cols1 = get_regex_cols(colnames(pheno_data),"pulse rate\\.",ignore.case=T)
additional_scores[['Pulse rate1']] = pheno_data[,pulse_rate_cols1]
pulse_rate_cols2 = get_regex_cols(colnames(pheno_data),"pulse rate,",ignore.case=T)
additional_scores[['Pulse rate2']] = pheno_data[,pulse_rate_cols2]
additional_scores[['Pulse rate']] = cbind(additional_scores[['Pulse rate1']],additional_scores[['Pulse rate2']])
hand_grip_cols = get_regex_cols(colnames(pheno_data),"Hand grip",ignore.case=T)
additional_scores[['hand grip']] = pheno_data[,hand_grip_cols]
rhr_v = additional_scores$`Pulse rate`
rhr_v = as.matrix(rhr_v)
mode(rhr_v) = 'numeric'
rhr_v = rowMeans(rhr_v,na.rm=T)
table(is.na(rhr_v))
rhr_v = rhr_v[!is.na(rhr_v)]
additional_scores[['Pulse rate']] = rhr_v

save(fitness_scores_matrix,additional_scores,accelerometry_scores,file="all_traits_for_gwas.RData")

###############################################
#################### End ######################
###############################################
###############################################

###############################################
###############################################
###### Normalize (quantile) and print #########
###############################################
###############################################
library(preprocessCore)

dir.create("uncorrected_fitness_scores/")
euro_pcs = read.delim("pca_results_v2_chrom1_euro.eigenvec")
euro_ids = as.character(euro_pcs$IID)
colnames(fitness_scores_matrix) = c("Recovery","HR_fitnes","Exercise_slopes","MaxWD","CompletionStatus")

# August 2017: print the scores as is
fitness_eu_ids = intersect(rownames(fitness_scores_matrix),euro_ids)
par(mfrow=c(4,2))
for(j in 1:ncol(fitness_scores_matrix)){
  currname = paste("uncorrected_fitness_scores/",colnames(fitness_scores_matrix)[j],"_euro.txt",sep='')
  m = cbind(fitness_eu_ids,fitness_eu_ids,fitness_scores_matrix[fitness_eu_ids,j])
  hist(as.numeric(m[,3]),main=paste(colnames(fitness_scores_matrix)[j],"uncorrected"))
  colnames(m) = c("FID","IID",colnames(fitness_scores_matrix)[j])
  write.table(m,file=currname,sep="\t",quote=F,row.names = F)
  m[,3] = gaus_norm(as.numeric(m[,3]))
  #plot(as.numeric(m[,3]),normalize.quantiles(as.matrix(as.numeric(m[,3])))) # sanity
  hist(as.numeric(m[,3]),main="normalized")
  currname = paste("uncorrected_fitness_scores/",colnames(fitness_scores_matrix)[j],"_qnorm_euro.txt",sep='')
  write.table(m,file=currname,sep="\t",quote=F,row.names = F)
}

rhr_v = additional_scores$`Pulse rate`
rhr_eu_samples = intersect(names(rhr_v),euro_ids)
currname = "uncorrected_fitness_scores/pulse_rate_euro.txt"
m = cbind(rhr_eu_samples,rhr_eu_samples,rhr_v[rhr_eu_samples])
colnames(m) = c("FID","IID","RHR")
write.table(m,file=currname,sep="\t",quote=F,row.names = F)
m[,3] = gaus_norm(as.numeric(m[,3]))
hist(as.numeric(m[,3]))
currname = "uncorrected_fitness_scores/pulse_rate_qnorm_euro.txt"
write.table(m,file=currname,sep="\t",quote=F,row.names = F)

###############################################
#################### End ######################
###############################################
###############################################

###############################################
###############################################
############## Filter features ################
###############################################
###############################################

# Define the sample set to work with: union of the above and rhr
subject_set = union(rownames(accelerometry_scores),rownames(fitness_scores_matrix))
subject_set = union(subject_set,names(additional_scores$`Pulse rate`))

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
print(table(percent_nas_col<1))
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

# Filter 6: exclude categorical features with too many classes (>50) 
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
table(features_to_exclude)
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
get_regex_cols(colnames(covariate_matrix),"batch",ignore.case=T)
get_regex_cols(colnames(covariate_matrix),"centr",ignore.case=T)
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

# Aug 2017: For now we use these only after correction
euro_pcs = read.delim("pca_results_v2_chrom1_euro.eigenvec")
euro_ids = as.character(euro_pcs$IID)

load("all_traits_for_gwas.RData")
rhr_v = additional_scores[["Pulse rate"]]
rm(additional_scores);gc()
colnames(fitness_scores_matrix) = c("Recovery","HR_fitnes","Exercise_slopes","MaxWD","CompletionStatus")
load("covariate_matrix.RData")
external_covs = read.delim('covariates.augmented.txt')
table(external_covs$PC1!=(-1000))
# exclude samples with NAs in their PCs: these do not have genotypes
external_covs = external_covs[external_covs$PC1!=(-1000),]
pcs_matrix = external_covs[,grepl("^PC",colnames(external_covs))]
rownames(pcs_matrix) = external_covs$IID
batch_info = as.character(external_covs$f.batch)
batch_info[batch_info==-1000] = NA
pcs_matrix[pcs_matrix==-1000] = NA
names(batch_info) = external_covs$IID
table(is.na(batch_info))
length(table(batch_info))
pcs_matrix = pcs_matrix[!apply(is.na(pcs_matrix),1,any),]
gc()
# # Compare Anna's covariates to the ones computed above
# sex1 = external_covs$Sex; names(sex1) = external_covs$IID
# sex2 = covariate_matrix[,"Sex"]; names(sex2) = rownames(covariate_matrix)
# sex_inters = intersect(names(sex1),names(sex2))
# table(sex1[sex_inters],sex2[sex_inters])
# Sex3: load the original pheno data
# load("biobank_collated_pheno_data.RData")
# Compare Anna's covariates to the ones computed above
# sex3 = pheno_data[,"Sex.0.0"]; names(sex3) = rownames(pheno_data)
# sex_inters = intersect(names(sex1),names(sex3))
# table(sex1[sex_inters],sex3[sex_inters])
# table(sex3[sex_inters],sex2[sex_inters])
# Compare the age
# age1 = external_covs$YearOfBirth;names(age1)=external_covs$IID
# age2 = covariate_matrix[,"Age when attended assessment centre"]
# names(age2) = rownames(covariate_matrix)
# sex_inters = sample(intersect(names(sex1),names(sex2)))[1:10000]
# plot(age1[sex_inters],age2[sex_inters])
# # Compare BMI
# bmi1 = external_covs$BMI;names(bmi1) = external_covs$IID
# bmi1[bmi1==-1000] = NA
# bmi2 = covariate_matrix[,"Body mass index (BMI)"];names(bmi2) = rownames(covariate_matrix)
# plot(bmi1[sex_inters],bmi2[sex_inters])
# merge the covariates and the pcs
names_inters = intersect(rownames(pcs_matrix),rownames(covariate_matrix))
covariate_matrix = cbind(covariate_matrix[names_inters,],pcs_matrix[names_inters,],batch_info[names_inters])
colnames(covariate_matrix)[ncol(covariate_matrix)] = "batch"
pc_names = colnames(pcs_matrix)
feature_is_numeric[pc_names]=T
feature_is_numeric["batch"] = F
covariate_matrix[,"batch"] = as.character(covariate_matrix[,"batch"])

# For printing a scores vector for gwas
print_scores_vector_for_gwas<-function(v,fname){
  nn = names(v)
  m = cbind(nn,nn,v)
  colnames(m) = c("FID","IID","RHR")
  write.table(m,file=fname,sep="\t",quote=F,row.names = F)
}

# For the fitness analysis only: get the category of the subjects
load("biobank_collated_pheno_data.RData")
category_matrix = pheno_data[rownames(fitness_scores_matrix),get_regex_cols(colnames(pheno_data),"category")]
cat_order = c("No category, ECG not to be done",
              "Category 4, at-rest measurement",
              "Category 3, cycle at constant level",
              "Category 2, cycle rising to 35% level",
              "Category 1, cycle rising to 50% level"
              )
cat_order = cat_order[length(cat_order):1]
get_disc_cat <- function(x,y){
  if(all(is.na(x))){return(NA)}
  return(which(is.element(y,set=x[!is.na(x)]))[1])
}
category_disc = apply(category_matrix,1,get_disc_cat,y=cat_order)
table(category_disc)
test_subjs = rownames(fitness_scores_matrix)[!is.na(fitness_scores_matrix[,2])]
table(category_disc[test_subjs])
plot(fitness_scores_matrix[test_subjs,2],category_disc[test_subjs])
rm(pheno_data); gc()

# Added on September 2017: set the sample set to europeans only
names_inters = intersect(rownames(covariate_matrix),euro_ids)
covariate_matrix = covariate_matrix[names_inters,]
fitness_scores_matrix = fitness_scores_matrix[intersect(rownames(fitness_scores_matrix),names_inters),]
rhr_v = rhr_v[intersect(names(rhr_v),names_inters)]

# Four possible GWAS runs:
# Simple vs Simpler analyses
# Outlier removal vs. qnorm
simple_covs = covariate_matrix[,c("Sex","Age when attended assessment centre","Standing height",pc_names,"batch")]
simple_covs[,"Sex"] = as.numeric(as.factor(simple_covs[,"Sex"]))-1
feature_is_numeric["Sex"]=T
simpler_covs = simple_covs[,c("Sex","Age when attended assessment centre",pc_names,"batch")]
fitness_vs_covs_lm_objects = list()
r2_table = c()
dir.create("gwas/simple_unnorm/");dir.create("gwas/simple_norm/")
dir.create("gwas/simpler_unnorm/");dir.create("gwas/simpler_norm/")
for (fitness_score_ind in 1:ncol(fitness_scores_matrix[,-1])){
  y = fitness_scores_matrix[,fitness_score_ind]
  y = y[!is.na(y)]
  y_stand = (y-mean(y,na.rm=T))/sd(y,na.rm=T)
  y_norm = gaus_norm(as.numeric(y));names(y_norm) = names(y)
  currname = colnames(fitness_scores_matrix)[fitness_score_ind]
  y = y[abs(y_stand)<6]
  par(mfrow=c(2,1));hist(y);hist(y_norm)
  lm1 = get_lm_residuals(y,simple_covs[names(y),],use_categorical=T,max_num_classes=150,feature_is_numeric=feature_is_numeric)
  lm2 = get_lm_residuals(y,simpler_covs[names(y),],use_categorical=T,max_num_classes=150,feature_is_numeric=feature_is_numeric)
  lm3 = get_lm_residuals(y_norm,simple_covs[names(y_norm),],use_categorical=T,max_num_classes=150,feature_is_numeric=feature_is_numeric)
  lm4 = get_lm_residuals(y_norm,simpler_covs[names(y_norm),],use_categorical=T,max_num_classes=150,feature_is_numeric=feature_is_numeric)
  res1 = lm1$lm_obj$residuals
  res2 = lm2$lm_obj$residuals
  res3 = lm3$lm_obj$residuals
  res4 = lm4$lm_obj$residuals
  print_scores_vector_for_gwas(res1,paste("gwas/simple_unnorm/",currname,".txt",sep=''))
  print_scores_vector_for_gwas(res2,paste("gwas/simpler_unnorm/",currname,".txt",sep=''))
  print_scores_vector_for_gwas(res3,paste("gwas/simple_norm/",currname,".txt",sep=''))
  print_scores_vector_for_gwas(res4,paste("gwas/simpler_norm/",currname,".txt",sep=''))
  r2s = c(
    get_lm_r2(lm1$lm_obj),get_lm_r2(lm2$lm_obj),get_lm_r2(lm3$lm_obj),get_lm_r2(lm4$lm_obj)
  )
  r2_table = rbind(r2_table,r2s)
  rownames(r2_table)[nrow(r2_table)] = currname
  #plot(y[names(lm1[[1]]$residuals)],lm1[[1]]$residuals,ylab="Residual",xlab="Fitness score")
  fitness_vs_covs_lm_objects[[currname]] = lm1
  save(fitness_vs_covs_lm_objects,file="fitness_analysis_fitness_vs_covs_lm_objects_simple.RData")
}
colnames(r2_table) = c("simple,unnorm","simpler,unnorm","simple,norm","simpler,norm")

# Add RHR to the simple norm
rhr_y = gaus_norm(rhr_v)
names(rhr_y) = names(rhr_v)
rhr_y = rhr_y[intersect(names(rhr_y),rownames(simpler_covs))]
currname = "RHR"
non_nas_cov_subjs = rownames(simple_covs)[apply(is.na(simple_covs),1,sum)==0]
rhr_y = rhr_y[intersect(names(rhr_y),non_nas_cov_subjs)]
rhr_lm1 = get_lm_residuals(rhr_y,simple_covs[names(rhr_y),],use_categorical=T,max_num_classes=150,feature_is_numeric=feature_is_numeric)
rhr_lm2 = get_lm_residuals(rhr_y,simpler_covs[names(rhr_y),],use_categorical=T,max_num_classes=150,feature_is_numeric=feature_is_numeric)
res3 = rhr_lm1$lm_obj$residuals
res3 = res3[intersect(euro_ids,names(res3))]
res4 = rhe_lm2$lm_obj$residuals
res4 = res4[intersect(euro_ids,names(res4))]
print_scores_vector_for_gwas(res3,paste("gwas/simple_norm/",currname,".txt",sep=''))
print_scores_vector_for_gwas(res4,paste("gwas/simpler_norm/",currname,".txt",sep=''))

# Plot for the report: Figure S2.2 
#load("fitness_analysis_fitness_vs_covs_lm_objects_simple.RData")
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

# new conservative analysis for fitness 
conservative_colnames = c("Sex","Age when attended assessment centre",
  "Standing height",pc_names,"batch","Smoking status",
  "Current employment status", "Alcohol drinker status",
  "Diabetes diagnosed by doctor", "Frequency of depressed mood in last 2 weeks",
  "Frequency of unenthusiasm / disinterest in last 2 weeks", 
  "Exposure to tobacco smoke at home", "Time spent using computer", "Time spent watching television (TV)"
)
conservative_covs = covariate_matrix[,conservative_colnames]
conservative_feature_type = feature_is_numeric[conservative_colnames]
conservative_covs[,"Sex"] = as.numeric(as.factor(conservative_covs[,"Sex"]))-1

fitness_vs_covs_lm_objects_conservative = list()
r2_cons = c()
dir.create("gwas/conservative_norm/")
for (fitness_score_ind in 1:ncol(fitness_scores_matrix[,-1])){
  y = fitness_scores_matrix[,fitness_score_ind]
  y = y[!is.na(y)]
  y_norm = gaus_norm(as.numeric(y));names(y_norm) = names(y)
  currname = colnames(fitness_scores_matrix)[fitness_score_ind]
  curr_covs = conservative_covs[names(y_norm),]
  curr_covs = curr_covs[,(colSums(is.na(curr_covs))/nrow(curr_covs)) < 0.2]
  lm1 = get_lm_residuals(y_norm,curr_covs,
                         use_categorical=T,max_num_classes=150,feature_is_numeric=conservative_feature_type)
  res1 = lm1$lm_obj$residuals
  res1 = res1[intersect(euro_ids,names(res1))]
  print_scores_vector_for_gwas(res1,paste("gwas/conservative_norm/",currname,".txt",sep=''))
  r2_cons[currname] = get_lm_r2(lm1$lm_obj)
  fitness_vs_covs_lm_objects_conservative[[currname]] = lm1
  save(fitness_vs_covs_lm_objects_conservative,file="fitness_analysis_fitness_vs_covs_lm_objects_conservative.RData")
}

# # accelerometry analysis
# accelerometry_scores_to_residuals = list()
# for (j in 1:ncol(accelerometry_scores)){
#   y = accelerometry_scores[,j]
#   names(y) = rownames(accelerometry_scores)
#   currname = colnames(accelerometry_scores)[j]
#   lm1 = get_lm_residuals(y,simple_covs,use_categorical=T,max_num_classes=20,feature_is_numeric=feature_is_numeric)
#   plot(y[names(lm1[[1]]$residuals)],lm1[[1]]$residuals,ylab="Residual",xlab="Fitness score")
#   accelerometry_scores_to_residuals[[currname]] = lm1[[1]]$residuals
#   save(accelerometry_scores_to_residuals,file="accelereometry_analysis_score_vs_covs_residuals_simple.RData")
# }

# Old conservative analysis
# fitness_vs_covs_lm_objects = list()
# for (fitness_score_ind in 1:ncol(fitness_scores_matrix)){
#   y = fitness_scores_matrix[,fitness_score_ind]
#   currname = colnames(fitness_scores_matrix)[fitness_score_ind]
#   lm1 = get_lm_residuals(y,covariate_matrix,use_categorical=T,max_num_classes=20,feature_is_numeric=feature_is_numeric)
#   plot(y[names(lm1[[1]]$residuals)],lm1[[1]]$residuals,ylab="Residual",xlab="Fitness score")
#   fitness_vs_covs_lm_objects[[currname]] = lm1
#   save(fitness_vs_covs_lm_objects,file="fitness_analysis_fitness_vs_covs_lm_objects.RData")
# }
# r2_scores = sapply(fitness_vs_covs_lm_objects,function(x)summary(x$lm)[["r.squared"]])

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

# Add the categorical accelerometry data
accelerometry_scores_discrete = read.delim(
  "accelerometry_phenotypes_anna/accelerometery_aggregate_phenotypes.categorical.filtered.txt",row.names = 1)
accelerometry_scores_discrete = accelerometry_scores_discrete[,-1]
colnames(accelerometry_scores_discrete)
# Simple
simple_covs = covariate_matrix[,c("Sex","Age when attended assessment centre","Body mass index (BMI)","Standing height")]
disc_accl_residual_scores = list()
for (j in 1:ncol(accelerometry_scores_discrete)){
  y = as.numeric(accelerometry_scores_discrete[,j])
  names(y) = rownames(accelerometry_scores_discrete)
  currname = colnames(accelerometry_scores_discrete)[j]
  lm1 = get_lm_residuals(y,simple_covs,use_categorical=T,max_num_classes=20,feature_is_numeric=feature_is_numeric,maxp=100000)
  plot(y[names(lm1[[1]]$residuals)],lm1[[1]]$residuals,ylab="Residual",xlab="Fitness score")
  disc_accl_residual_scores[[currname]] = lm1[[1]]$residuals
  save(disc_accl_residual_scores,file="disc_accl_residual_scores_simple.RData")
}
# Conservative
disc_accl_residual_scores = list()
for (j in 1:ncol(accelerometry_scores_discrete)){
  y = as.numeric(accelerometry_scores_discrete[,j])
  names(y) = rownames(accelerometry_scores_discrete)
  currname = colnames(accelerometry_scores_discrete)[j]
  lm1 = get_lm_residuals(y,covariate_matrix,use_categorical=T,max_num_classes=20,feature_is_numeric=feature_is_numeric)
  plot(y[names(lm1[[1]]$residuals)],lm1[[1]]$residuals,ylab="Residual",xlab="Fitness score")
  disc_accl_residual_scores[[currname]] = lm1[[1]]$residuals
  save(disc_accl_residual_scores,file="disc_accl_residual_scores_conservative.RData")
}

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
analyze_associations_between_scores_and_covariates<-function(scores,subject_set,cov_mat,disc_cov_mat,cov2category){
  summary_table = c()
  fv = scores[subject_set]
  fv_disc = cut_by_quantiles(fv,nbreaks = cut_size)
  fv_is_na = is.na(fv)
  for(j in 1:ncol(cov_mat)){
    covv = cov_mat[subject_set,j]
    covv_disc = Mc[subject_set,j]
    covv_na = is.na(covv)
    cov_name = colnames(cov_mat)[j]
    cov_category = cov2category[cov_name]
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
    cov_summary_scores = c(disc_MI,disc_chisq_pvalue,fv_na_mi,fv_na_chisq_p,
                           sp_rho,sp_rho_p,na_vs_na_MI,na_vs_na_chisq_pvalue)
    summary_table = rbind(summary_table,c(cov_name,cov_category,cov_summary_scores))
  }
  colnames(summary_table) = c("Feature","Category","MI-discretized","ChisqP-discretized",
                              "MI-NA fitness","ChisqP-NA fitness","Spearman rho","Spearman rho p","MI-NAs", "ChisqP-NAs")
  return(summary_table)
}
# scores = fitness_scores_matrix[,1]
# names(scores) = rownames(fitness_scores_matrix)
# tmp = analyze_associations_between_scores_and_covariates(scores,curr_subject_set,covariate_matrix,Mc,feature_category_data)
# tmp[1:3,] == covariance_correlation_summary_tables[[1]][1:3,]

# Analyze the fitness scores
covariance_correlation_summary_tables = list()
curr_subject_set = intersect(rownames(covariate_matrix),rownames(fitness_scores_matrix))
for(fit_ind in 1:ncol(fitness_scores_matrix)){
  scores = fitness_scores_matrix[,fit_ind]
  names(scores) = rownames(fitness_scores_matrix)
  summary_table = analyze_associations_between_scores_and_covariates(scores,curr_subject_set,covariate_matrix,Mc,feature_category_data)
  ord = order(as.numeric(summary_table[,3]),decreasing=T)
  print(summary_table[ord[1:5],1:4])
  covariance_correlation_summary_tables[[colnames(fitness_scores_matrix)[fit_ind]]] = summary_table
}
names(covariance_correlation_summary_tables) = colnames(fitness_scores_matrix)
save(covariance_correlation_summary_tables,file="fitness_analysis_covariance_correlation_summary_tables.RData")

summary_table_for_the_stext2 = c()
out_file = "SText2_STable2.1.txt"
table_cols = c(1,2,3,4,7)
for(j in 1:length(covariance_correlation_summary_tables)){
  nn = names(covariance_correlation_summary_tables)[j]
  summary_table = covariance_correlation_summary_tables[[nn]]
  is_max_in_cat = rep(F,nrow(summary_table))
  for(i in 1:nrow(summary_table)){
    curr_cat = summary_table[i,2]
    curr_max = max(as.numeric(summary_table[summary_table[,2]==curr_cat,3]),na.rm=T)
    if(as.numeric(summary_table[i,3])==curr_max){
      is_max_in_cat[i]=T
    }
  }
  summary_table = summary_table[is_max_in_cat,]
  ord = order(as.numeric(summary_table[,3]),decreasing=T)
  summary_table = summary_table[ord,table_cols]
  summary_table = cbind(rep(nn,nrow(summary_table)),summary_table)
  summary_table_for_the_stext2 = rbind(summary_table_for_the_stext2,summary_table)
  write.table(summary_table,file=out_file,sep="\t",quote=F,row.names = F,col.names = j==1,append = j>1)
}

# Analyze the accelerometry
acc_mat = accelerometry_scores
rows_with_nas = apply(is.na(acc_mat),1,any)
total_na_percent_colm = colSums(is.na(acc_mat))/nrow(acc_mat)
table(total_na_percent_colm)
acc_mat = acc_mat[,total_na_percent_colm<0.8]
total_na_percent_row = rowSums(is.na(acc_mat))/ncol(acc_mat)
acc_mat = as.matrix(acc_mat[total_na_percent_row<0.5,])
acc_mat = impute.knn(acc_mat,rowmax = 0.8)$data
acc_pca = prcomp(acc_mat,retx=T)
plot(acc_pca)
acc_pc1 = acc_pca$x[,1]
curr_subject_set = intersect(rownames(covariate_matrix),rownames(acc_mat))
acc_vs_covariates = analyze_associations_between_scores_and_covariates(acc_pc1,curr_subject_set,covariate_matrix,Mc,feature_category_data)
ord = order(as.numeric(acc_vs_covariates[,3]),decreasing=T)
print(acc_vs_covariates[ord[1:10],1:4])
write.table(acc_vs_covariates,file="Accelerometry_data_PC1_vs_covariates.txt",quote=F,row.names=F,sep="\t")

# Look at NAs vs. non NAs scores
names(covariance_correlation_summary_tables)
ind=2
currname = names(covariance_correlation_summary_tables)[ind]
summary_table = covariance_correlation_summary_tables[[ind]]
colnames(summary_table)

###############################################
###############################################
#################### End ######################
###############################################
###############################################

###############################################
###############################################
############# Tests and Figures ###############
###############################################
###############################################

# Figure S2.1
x1 = -log(pmax(1e-200,as.numeric(summary_table[,"ChisqP-discretized"])),10)
y1 = -log(pmax(1e-200,as.numeric(summary_table[,"ChisqP-NA fitness"])),10)
par(mfrow=c(1,2),mar=c(5,5,5,5))
plot(x=x1,y=y1,
     main="Predicted HR: p-value",xlab="-log p scores",ylab="-log p NAs",pch=4,lwd=1.5);abline(0,1)
summary_table[which(y1>150 & x1<10),1]
summary_table[which(y1==200 & x1==200),1]
plot(x=summary_table[,"MI-discretized"],y=summary_table[,"MI-NA fitness"],
     main="Predicted HR: mutual information",xlab="MI vs. scores",ylab="MI vs. NAs",pch=4,lwd=1.5,xlim=c(0,0.2),ylim=c(0,0.2));abline(0,1)

# Compare our results to VERSION 1
load("UKBB_phenotypic_data_for_GWAS.RData")
x1 = fitness_scores_matrix[,2]
colnames(fitness_scores)
x2 = fitness_scores[,2]
curr_inter = intersect(names(x1),names(x2))
get_pairwise_corrs(cbind(x1[curr_inter],x2[curr_inter]),method="spearman")

