# In this script we combine all phenotypes and put them in matrices to be analyzed with plink

try({setwd("/Users/davidhsu/Documents/ukbb")})
try({setwd("/Users/david/Desktop/ukbb/")})
try({setwd('/scratch/PI/euan/projects/ukbb/da_dh/')})

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
save(final_mat,file="physical_activity_scores_for_GWAS.RData")

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
save(final_mat,file="additional_scores_for_GWAS.RData")

###############################################
###############################################
#################### End ######################
###############################################
###############################################