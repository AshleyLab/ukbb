# In this script we combine all phenotypes and put them in matrices to be analyzed with plink

try({setwd("/Users/davidhsu/Documents/ukbb")})
try({setwd("/Users/david/Desktop/ukbb/")})
try({setwd('/scratch/PI/euan/projects/ukbb/da_dh/')})
source("auxiliary_functions.R")
library(corrplot)

###############################################
###############################################
############# Helper functions ################
###############################################
###############################################

get_score_from_object<-function(x){
  if(class(x)!="numeric"){
    return(x[[1]]$residuals)
  }
  return(x)
}

# get_binary_score_by_NAs<-function(x){
#   v = rep(0,length(x))
#   names(v) = names(x)
#   v[!is.na(x)] = 1
#   return(v)
# }

merge_scores_list_into_matrix<-function(l){
  allnames = unique(unlist(sapply(l,names)))
  m = matrix(NA,nrow=length(allnames),ncol=length(l))
  colnames(m) = names(l)
  rownames(m) = allnames
  for(j in 1:ncol(m)){
    x = l[[j]]
    m[names(x),j]=x
  }
  return(m)
}

###############################################
###############################################
#################### End ######################
###############################################
###############################################

###############################################
###############################################
########## Physical fitness scores ############
###############################################
###############################################

# Conservative
file = "fitness_analysis_fitness_vs_covs_lm_objects.RData"
obj = get(load(file))
scores = lapply(obj,get_score_from_object)
scores_mat_consv = merge_scores_list_into_matrix(scores)
colnames(scores_mat_consv) = paste(colnames(scores_mat_consv),"_conservative",sep='')

# Simple
file = "fitness_analysis_fitness_vs_covs_lm_objects_simple.RData"
obj = get(load(file))
scores = lapply(obj,get_score_from_object)
scores_mat_simp = merge_scores_list_into_matrix(scores)
colnames(scores_mat_simp) = paste(colnames(scores_mat_simp),"_simple",sep='')

final_mat = cbind(scores_mat_simp,scores_mat_consv)
print(dim(final_mat))
corrs = get_pairwise_corrs(final_mat)
corrplot(corrs,order='hclust')
save(final_mat,file="physical_fitness_scores_for_GWAS.RData")

###############################################
###############################################
#################### End ######################
###############################################
###############################################

###############################################
###############################################
########## Physical activity scores ###########
###############################################
###############################################

# Conservative
file = "accelereometry_analysis_score_vs_covs_residuals_conservative.RData"
obj = get(load(file))
scores = lapply(obj,get_score_from_object)
scores_mat_consv = merge_scores_list_into_matrix(scores)
colnames(scores_mat_consv) = paste(colnames(scores_mat_consv),"_conservative",sep='')

# Simple
file = "accelereometry_analysis_score_vs_covs_residuals_simple.RData"
obj = get(load(file))
scores = lapply(obj,get_score_from_object)
scores_mat_simp = merge_scores_list_into_matrix(scores)
colnames(scores_mat_simp) = paste(colnames(scores_mat_simp),"_simple",sep='')

final_mat = cbind(scores_mat_simp,scores_mat_consv)
print(dim(final_mat))
corrs = get_pairwise_corrs(final_mat[,1:20])
corrplot(corrs,order='hclust')
save(final_mat,file="physical_activity_scores_for_GWAS.RData")

# The categorical data
# Conservative
file = "disc_accl_residual_scores_conservative.RData"
obj = get(load(file))
scores = lapply(obj,get_score_from_object)
scores_mat_consv = merge_scores_list_into_matrix(scores)
colnames(scores_mat_consv) = paste(colnames(scores_mat_consv),"_conservative",sep='')
# Simple
file = "disc_accl_residual_scores_simple.RData"
obj = get(load(file))
scores = lapply(obj,get_score_from_object)
scores_mat_simp = merge_scores_list_into_matrix(scores)
colnames(scores_mat_simp) = paste(colnames(scores_mat_simp),"_simple",sep='')
final_mat = cbind(scores_mat_simp,scores_mat_consv)
print(dim(final_mat))
corrs = get_pairwise_corrs(final_mat[,1:20])
corrplot(corrs,order='hclust')
save(final_mat,file="discrete_physical_activity_scores_for_GWAS.RData")

###############################################
###############################################
#################### End ######################
###############################################
###############################################

###############################################
###############################################
############## Additional scores ##############
###############################################
###############################################

# Conservative
file = "additional_scores_vs_covs_lm_objects_vs_covs_lm_objects_conservative.RData"
obj = get(load(file))
scores = lapply(obj,get_score_from_object)
scores_mat_consv = merge_scores_list_into_matrix(scores)
colnames(scores_mat_consv) = paste(colnames(scores_mat_consv),"_conservative",sep='')

# Simple
file = "additional_scores_vs_covs_lm_objects_vs_covs_lm_objects_simple.RData"
obj = get(load(file))
scores = lapply(obj,get_score_from_object)
scores_mat_simp = merge_scores_list_into_matrix(scores)
colnames(scores_mat_simp) = paste(colnames(scores_mat_simp),"_simple",sep='')

final_mat = cbind(scores_mat_simp,scores_mat_consv)
print(dim(final_mat))
corrs = get_pairwise_corrs(final_mat)
corrplot(corrs,order='hclust')
save(final_mat,file="additional_scores_for_GWAS.RData")

###############################################
###############################################
#################### End ######################
###############################################
###############################################

###############################################
###############################################
################ Some plots ###################
###############################################
###############################################

library(corrplot)
corrs = get_pairwise_corrs(final_mat)
corrplot(corrs,order='hclust')

# Merge matrices
l = list()
l[["additional"]] = get(load("additional_scores_for_GWAS.RData"))
l[["fitness"]] = get(load("physical_fitness_scores_for_GWAS.RData"))
all_pheno_mat = merge_scores_matriceslist_into_matrix(l)
dim(all_pheno_mat)
save(all_pheno_mat,file="all_pheno_mat_for_GWAS.RData")

# Assumes that column names across all matrices are unique
merge_scores_matriceslist_into_matrix<-function(l){
  rows = unique(unlist(sapply(l,rownames)))
  cols = unlist(sapply(l,colnames))
  m = matrix(NA,nrow=length(rows),ncol=length(cols))
  colnames(m) = cols
  rownames(m) = rows
  for(j in 1:length(l)){
    x = l[[j]]
    m[rownames(x),colnames(x)]=x
  }
  return(m)
}

# Sanity checks
table(is.na(all_pheno_mat[,15]))
table(is.na(l[[2]][,5]))
x1 = all_pheno_mat[,15]
x2 = l[[2]][,5]
inter = intersect(names(x1),names(x2))
all(x1[inter]==x2[inter],na.rm=T)
all(is.na(x1[inter])==is.na(x2[inter]))
plot(x1[inter],x2[inter])

###############################################
###############################################
#################### End ######################
###############################################
###############################################

###############################################
###############################################
############## Correct PCs ####################
###############################################
###############################################
source("auxiliary_functions.R")
scores_matrix = get(load("physical_fitness_scores_for_GWAS.RData"))

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
# Correct for the PCs
for(nn in colnames(scores_matrix)){
  curr_res = scores_matrix[,nn]
  curr_res = curr_res[gwas_data_samples]
  names(curr_res) = gwas_data_samples
  NA_samples = gwas_data_samples[is.na(curr_res)]
  gwas_data_residuals = cbind(gwas_data_residuals,curr_res)
}
colnames(gwas_data_residuals) = paste("Residuals_",colnames(scores_matrix),sep="")
colnames(gwas_data_residuals) = gsub(colnames(gwas_data_residuals),pattern=" ",replace="_")
gwas_data = cbind(gwas_data_residuals,as.matrix(gwas_data_pcs))
corrs = get_pairwise_corrs(gwas_data)
save(gwas_data,file="July19_2017_gwas_data_table.RData")
write.table(gwas_data,file = "July19_2017_fitness_scores_gwas_data_table.txt",sep="\t",quote=F)
library(corrplot)
corrplot(corrs,order='hclust')





