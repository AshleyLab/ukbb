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
x = as.character(x)
x[x=="Category 1, cycle rising to 50% level"] = "\nCategory 1\ncycle to 50% level"
x[x=="Category 2, cycle rising to 35% level"] = "\nCategory 2\ncycle to 35% level"
x = factor(x)
boxplot(y~x,ylab="Recovery (percentage)", main = "Recovery",notch=T,xaxt = "n")
axis(side = 1,at = 1:2,levels(x),tick = FALSE)

# Plot the scores vs. the categories: Fitness
y = fitness_scores_matrix[,2]
y = y[!is.na(y)]
x = subj2category[names(y)]
x = x[grepl("Category 1|2",x)]
y = y[names(x)]
y = y[abs((y-mean(y))/sd(y))<6]
x = x[names(y)]
x = as.character(x)
x[x=="Category 1, cycle rising to 50% level"] = "\nCategory 1\ncycle to 50% level"
x[x=="Category 2, cycle rising to 35% level"] = "\nCategory 2\ncycle to 35% level"
x = factor(x)
boxplot(y~x,ylab="Exercise HR",main="Exercise HR", notch=T,xaxt = "n")
axis(side = 1,at = 1:2,levels(x),tick = FALSE) 

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
library(hexbin)
plot_two_specific_scores<-function(x1,x2,...){
  x1 = x1[!is.na(x1)]
  x2 = x2[!is.na(x2)]
  x1 = x1[abs((x1-mean(x1))/sd(x1))<6]
  x2 = x2[abs((x2-mean(x2))/sd(x2))<6]
  nns = intersect(names(x1),names(x2))
  x1 = x1[nns];x2 = x2[nns]
  r2 = cor(x1,x2)^2
  sp = cor(x1,x2,method="spearman")
  gplot.hexbin(hexbin(x1,x2,xbins=50,...),
               legend=F,main=paste("R^2=",format(r2,digits = 3),", Spearman rho=",format(sp,digits=3),sep=""))
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
# names(additional_scores)
# additional_scores[[1]] = additional_scores[[1]][,105:106]
# additional_scores[[3]] = additional_scores[[3]][,c(6,19)]
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
additional_scores_cov_analysis = list()
additional_scores_cov_analysis[["ExerciseHR"]] = analyze_associations_between_scores_and_covariates(curr_scores,
      additional_scores_mat,additional_scores_mat_disc1,feature_subcategory_data,add_scores_f2isnum)
curr_scores = fitness_scores_matrix[,1]
curr_scores = curr_scores[!is.na(curr_scores)]
additional_scores_cov_analysis[["Recovery"]] = analyze_associations_between_scores_and_covariates(curr_scores,
      additional_scores_mat,additional_scores_mat_disc1,feature_subcategory_data,add_scores_f2isnum)


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

save(cov_scores_recovery,cov_scores_hr,additional_scores_cov_analysis,file="fitness_scores_mi_rho_vs_other_scores.RData")

load("fitness_scores_mi_rho_vs_other_scores.RData")
# Plot the associations
col_inds = c(1,2,3,5)
# recovery
scores1 = cov_scores_recovery[,col_inds]
scores2 = additional_scores_cov_analysis$Recovery[,col_inds]
scores = rbind(scores1,scores2)
ord = order(abs(as.numeric(scores[,3])),decreasing=T,na.last = T)
scores = scores[ord,]
rownames(scores) = scores[,1]
scores = data.frame(scores)
scores[,3] = as.numeric(as.character(scores[,3]))
scores[,4] = as.numeric(as.character(scores[,4]))
plot(scores$MI.discretized,scores$Spearman.rho)
sub_score_names = tapply(scores[,1],scores[,2],function(x)names(x)[1])
sub_scores = data.frame(scores[sub_score_names,])
sub_scores = sub_scores[sub_scores$MI.discretized>0.001,]
dim(sub_scores)
rec_sub_scores = sub_scores

cov_scores_recovery[grepl("hae",cov_scores_recovery[,1]),]


# Exercise
scores1 = cov_scores_hr[,col_inds]
scores2 = additional_scores_cov_analysis$ExerciseHR[,col_inds]
scores = rbind(scores1,scores2)
ord = order(abs(as.numeric(scores[,3])),decreasing=T,na.last = T)
scores = scores[ord,]
rownames(scores) = scores[,1]
scores = data.frame(scores)
scores[,3] = as.numeric(as.character(scores[,3]))
scores[,4] = as.numeric(as.character(scores[,4]))
plot(scores$MI.discretized,scores$Spearman.rho)
sub_score_names = tapply(scores[,1],scores[,2],function(x)names(x)[1])
sub_scores = data.frame(scores[sub_score_names,])
sub_scores = sub_scores[sub_scores$MI.discretized>0.001,]
dim(sub_scores)
ex_sub_scores = sub_scores

rec_sub_scores = cbind(rep("Recovery",nrow(rec_sub_scores)),rec_sub_scores)
ex_sub_scores = cbind(rep("ExerciseHR",nrow(ex_sub_scores)),ex_sub_scores)
all_sub_scores = rbind(as.matrix(rec_sub_scores),as.matrix(ex_sub_scores))
rownames(all_sub_scores) = NULL
colnames(all_sub_scores)[1] = "Pheno"
ord = order(as.numeric(all_sub_scores[,4]),decreasing = T)
all_sub_scores = all_sub_scores[ord,]
write.table(all_sub_scores,file="fitness_analysis_covariance_correlation_table.txt",sep="\t",quote=F)

# Look at some interesting correlations
# 1. Hemoglobin
par(mfrow=c(1,2))
x1 = covariate_matrix[,"Haemoglobin concentration"];names(x1) = rownames(covariate_matrix)
x2 = rhr_v
nns = intersect(names(x1),names(x2))
x1 = x1[nns];x2=x2[nns]
plot_two_specific_scores(x1,x2,xlab="Hemoglobin concentration",ylab="Resting heart rate")
x2 = fitness_scores_matrix[,2]
nns = intersect(names(x1),names(x2))
x1 = x1[nns];x2=x2[nns]
plot_two_specific_scores(x1,x2,xlab="Hemoglobin concentration",ylab="Exercise HR")

scores_for_pairwise_analysis = list()
x1 = covariate_matrix[,"Haemoglobin concentration"];names(x1) = rownames(covariate_matrix)
scores_for_pairwise_analysis[["Hemoglobin"]] = x1
scores_for_pairwise_analysis[["ExerciseHR"]] = fitness_scores_matrix[,2]
scores_for_pairwise_analysis[["Recovery"]] = fitness_scores_matrix[,1]
scores_for_pairwise_analysis[["RHR"]] = additional_scores$`Pulse rate`
scores_for_pairwise_analysis[["HandGrip"]] = additional_scores_mat[,'hand grip']
x1 = covariate_matrix[,"Body mass index (BMI)"];names(x1) = rownames(covariate_matrix)
scores_for_pairwise_analysis[["BMI"]] = x1
x1 = covariate_matrix[,"Age when attended assessment centre"];names(x1) = rownames(covariate_matrix)
scores_for_pairwise_analysis[["Age"]] = x1
x1 = covariate_matrix[,"Body fat percentage"];names(x1) = rownames(covariate_matrix)
scores_for_pairwise_analysis[["Fat percentage"]] = x1
x1 = covariate_matrix[,"Standing height"];names(x1) = rownames(covariate_matrix)
scores_for_pairwise_analysis[["Height"]] = x1
x1 = covariate_matrix[,"Number of days/week of vigorous physical activity 10+ minutes"];names(x1) = rownames(covariate_matrix)
scores_for_pairwise_analysis[["Self reported"]] = x1
corrs = get_pairwise_corrs(scores_for_pairwise_analysis)
corrplot(corrs)

"Diastolic blood pressure, automated reading"
"Systolic blood pressure, automated reading"
"Alcohol intake frequency"
"Time spent using computer"
"Number of days/week of vigorous physical activity 10+ minutes" 

# transform corrs to table for visualization in cytoscape
mat = c()
for(i in 2:ncol(corrs)){
  for(j in 1:(i-1)){
    n1 = colnames(corrs)[i]
    n2 = colnames(corrs)[j]
    mat = rbind(mat,c(n1,n2,corrs[i,j]))
  }
}
write.table(mat,file="fitness_analysis_selected_covariance_correlation_table.txt",
            sep="\t",quote=F,row.names = F,col.names = F)

# Conditional independence tests
library(bnlearn)
# d must have x,y, and z
run_discrete_ci_test<-function(x,y,z,cutsize=10,contin=T,...){
  nns = intersect(names(x),names(y))
  nns = intersect(nns,names(z))
  d = data.frame(x=x[nns],y=y[nns],z=z[nns])
  d = d[!apply(is.na(d),1,any),]
  if(contin){return(ci.test(d,...))}
  for(nn in names(d)){
    if(is.numeric(d[[nn]])){
      d[[nn]] = factor(cut(d[[nn]],breaks=cutsize))
    }
  }
  return(ci.test(d,...))
}
run_discrete_ci_test(scores_for_pairwise_analysis[[2]],
                     scores_for_pairwise_analysis[[1]],
                     scores_for_pairwise_analysis[[5]],test='cor')


