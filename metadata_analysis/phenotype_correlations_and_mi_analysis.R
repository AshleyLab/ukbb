# The following script loads the phenotypic data
# We start by running some QA tests in order to characterize the experiment that was
# taken for each individual
# On sherlock use: module load R

try({setwd("/Users/davidhsu/Documents/ukbb")})
try({setwd("/Users/david/Desktop/ukbb/")})
try({setwd('/scratch/PI/euan/projects/ukbb/da_dh/')})
source("auxiliary_functions.R")

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

# Fitness scores:
load("all_traits_for_gwas.RData")
load("covariate_matrix.RData")

# Inspect the subject categories, are there discrepancies?
category_matrix = pheno_data[,get_regex_cols(colnames(pheno_data),"category")]
is_same_category<-function(x){
  if(any(is.na(x))){return(T)}
  return (x[1]==x[2])
}
subj2category = category_matrix[,1]
names(subj2category) = rownames(category_matrix) 
subj2category[is.na(subj2category)] = category_matrix[is.na(subj2category),2]  
table(subj2category)

# Plot the scores vs. the categories: Recovery
y = fitness_scores_matrix[,1]
y = y[!is.na(y)]
x = subj2category[names(y)]
x = x[grepl("Category 1|2",x)]
y = y[names(x)]
y = y[abs((y-mean(y))/sd(y))<6]
x = x[names(y)]
x = factor(as.character(x))
boxplot(y~x,ylab="Recovery (percentage)", main = "Recovery")
x = as.numeric(x)
wilcox.test(y[x==1],y[x==2])$p.value  


# Plot the scores vs. the categories: Fitness
y = fitness_scores_matrix[,2]
y = y[!is.na(y)]
x = subj2category[names(y)]
x = x[grepl("Category 1|2",x)]
y = y[names(x)]
y = y[abs((y-mean(y))/sd(y))<6]
x = x[names(y)]
x = factor(as.character(x))
boxplot(y~x,ylab="Exercise HR",main="Exercise HR")  
x = as.numeric(x)
wilcox.test(y[x==1],y[x==2])$p.value

# Aux functions
cut_by_quantiles<-function(x,nbreaks=5){
  if(length(unique(x))<=nbreaks){return(x)}
  probs = seq(0,1,length.out = nbreaks+1)
  qs = unique(quantile(x,probs = probs,na.rm=T))
  xv = cut(x,breaks=qs)
  return(xv)
}

library(entropy)
analyze_associations_between_scores_and_covariates<-function(scores,cov_mat,disc_cov_mat,cov2category,feature_is_numeric){
  summary_table = c()
  subject_set = names(scores)
  fv = scores
  fv_disc = cut_by_quantiles(fv,nbreaks = cut_size)
  fv_is_na = is.na(fv)
  for(j in 1:ncol(cov_mat)){
    covv = cov_mat[subject_set,j]
    covv_disc = disc_cov_mat[subject_set,j]
    covv_na = is.na(covv)
    cov_name = colnames(cov_mat)[j]
    cov_category = cov2category[cov_name]
    # correlations
    non_na_inds = !fv_is_na & !is.na(covv)
    sp_rho = NA ; sp_rho_p = NA; r2=NA
    if(sum(non_na_inds)>1000 && feature_is_numeric[cov_name]){
      x1 = as.numeric(covv[non_na_inds])
      x2 = fv[non_na_inds]
      sp_rho = cor(x1,x2,method='spearman')
      sp_rho_p = cor.test(x1,x2,method='spearman')$p.value
      r2 = summary(lm(x2~x1))[["r.squared"]]
    }
    # discrete vs discrete
    disc_chisq_pvalue = NA;disc_MI = NA
    if(length(unique(covv_disc))<100){
      tab = table(as.character(covv_disc),fv_disc)
      disc_chisq_pvalue = chisq.test(tab)$p.value
      disc_MI = mi.empirical(tab)
    }
    cov_summary_scores = c(disc_MI,disc_chisq_pvalue,sp_rho,sp_rho_p)
    summary_table = rbind(summary_table,c(cov_name,cov_category,cov_summary_scores,r2))
  }
  colnames(summary_table) = c("Feature","Category","MI-discretized","ChisqP-discretized"
                              ,"Spearman rho","Spearman rho p","R2")
  return(summary_table)
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

# Look at correlations and MI
names(additional_scores)
additional_scores[[1]] = additional_scores[[1]][,105:106]
additional_scores[[3]] = additional_scores[[3]][,c(6,19)]
sapply(additional_scores,colnames)
add_scores_subjs = unique(unlist(sapply(additional_scores,names)))
add_scores_subjs = union(add_scores_subjs,unique(unlist(sapply(additional_scores,rownames))))
additional_scores_mat = matrix(NA,ncol=5,nrow=length(add_scores_subjs))
rownames(additional_scores_mat) = add_scores_subjs
colnames(additional_scores_mat) = names(additional_scores)[c(1:3,6:7)]
specific_names = c()
for(cc in colnames(additional_scores_mat)){
  x = additional_scores[[cc]]
  if(!is.null(dim(x))){
    colnames(x)  = sapply(colnames(x),function(x)strsplit(x,split='\\.')[[1]][1])
    num_vals = apply(!is.na(x),2,sum)
    spec_name = names(which(num_vals==max(num_vals)))[1]
    x = x[,spec_name];names(x) = rownames(additional_scores[[cc]])
    specific_names = c(specific_names,spec_name)
  }
  else{
    specific_names = c(specific_names,spec_name)
  }
  additional_scores_mat[names(x),cc] = as.numeric(x)
}
par(mfrow=c(2,3))
apply(additional_scores_mat,2,hist)
additional_scores_mat_disc1 = apply(additional_scores_mat,2,cut_by_quantiles)
rownames(additional_scores_mat_disc1) = rownames(additional_scores_mat)
additional_scores_mat_disc2 = apply(additional_scores_mat,2,cut,breaks=5)
rownames(additional_scores_mat_disc2) = rownames(additional_scores_mat)
apply(additional_scores_mat_disc2,2,table)
apply(additional_scores_mat_disc1,2,table)

add_scores_f2isnum = rep(T,ncol(additional_scores_mat_disc1))
names(add_scores_f2isnum) = colnames(additional_scores_mat_disc1)
curr_scores = fitness_scores_matrix[,2]
curr_scores = curr_scores[!is.na(curr_scores)]
analyze_associations_between_scores_and_covariates(curr_scores,
      additional_scores_mat,additional_scores_mat_disc1,feature_subcategory_data,add_scores_f2isnum)

plot_two_specific_scores<-function(x1,x2,...){
  x1 = x1[!is.na(x1)]
  x2 = x2[!is.na(x2)]
  x1 = x1[abs((x1-mean(x1))/sd(x1))<6]
  x2 = x2[abs((x2-mean(x2))/sd(x2))<6]
  nns = intersect(names(x1),names(x2))
  x1 = x1[nns];x2 = x2[nns]
  r2 = cor(x1,x2)^2
  sp = cor(x1,x2,method="spearman")
  plot(x1,x2,main=paste("R^2=",format(r2,digits = 3),", Spearman rho=",format(sp,digits=3),sep=""),...)
}
par(mfrow=c(2,2))
curr_scores = fitness_scores_matrix[,2]
plot_two_specific_scores(curr_scores,additional_scores_mat[,1],ylab="Hand grip",xlab="Exercise HR")
plot_two_specific_scores(curr_scores,additional_scores_mat[,5],ylab="Resting HR",xlab="Exercise HR")
curr_scores = fitness_scores_matrix[,1]
plot_two_specific_scores(curr_scores,additional_scores_mat[,1],ylab="Hand grip",xlab="Recovery")
plot_two_specific_scores(curr_scores,additional_scores_mat[,5],ylab="Resting HR",xlab="Recovery")
par(mfrow=c(1,1))
plot_two_specific_scores(fitness_scores_matrix[,1],fitness_scores_matrix[,2],ylab="Exercise HR",xlab="Recovery")

curr_scores = fitness_scores_matrix[,1]
curr_scores = curr_scores[!is.na(curr_scores)]
cov_scores_recovery = analyze_associations_between_scores_and_covariates(curr_scores,
    covariate_matrix,Mc,feature_subcategory_data,feature_is_numeric)
curr_scores = fitness_scores_matrix[,2]
curr_scores = curr_scores[!is.na(curr_scores)]
cov_scores_hr = analyze_associations_between_scores_and_covariates(curr_scores,
    covariate_matrix,Mc,feature_subcategory_data,feature_is_numeric)

ord = order(as.numeric(cov_scores_recovery[,3]),decreasing=T,na.last = T)
cov_scores_recovery[ord,][1:5,]
ord = order(as.numeric(cov_scores_hr[,3]),decreasing=T,na.last = T)
cov_scores_hr[ord,][1:50,]
















