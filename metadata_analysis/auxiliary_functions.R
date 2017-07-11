get_regex_cols<-function(cols,reg,...){return (cols[grepl(cols,pattern=reg,...)])}
get_list_by_values<-function(x,values,min_size=1,values_to_ignore = NULL){
  l = list()
  val_set = sort(unique(values))
  if(!is.null(values_to_ignore)){
    val_set = setdiff(val_set,values_to_ignore)
  }
  for(v in val_set){
    inds = values==v
    if(sum(inds)<min_size){next}
    l[[as.character(v)]] = x[inds]
  }
  return(l)
}
#Return the names of the columns of a particular matrix that match the desired pattern.
get_regex_cols<-function(cols,reg,...){return (cols[grepl(cols,pattern=reg,...)])}
get_last_index_that_fits_value<-function(x,val){
  val_inds = which(x==val)
  if(length(val_inds)==0){return(NA)}
  return (val_inds[length(val_inds)])
}
# This function returns the first index within val that matches x.
get_first_index_that_fits_value<-function(x,val){
  val_inds = which(x==val)
  if(length(val_inds)==0){return(NA)}
  return(val_inds[1])
}
run_lm<-function(y,x,...){
  obj = lm(y~.,data=data.frame(y,x))
  return(obj)
}

#install.packages('hexbin',lib = '/home/davidama/R/x86_64-pc-linux-gnu-library/3.2/')
try({library(hexbin,lib.loc='/home/davidama/R/x86_64-pc-linux-gnu-library/3.2/')})
try({library(hexbin)})
# A set of functions for dealing with vectors with different name set
pairwise_plot<-function(x,y,func = joint.density.plot,sample_set = NULL,...){
  intr = intersect(names(x),names(y))
  if(!is.null(sample_set)){intr=intersect(sample_set,intr)}
  func(x[intr],y[intr],...)
}
# R^2 based analysis
get_lm_r2<-function(x)unlist(unname(summary(x)["r.squared"]))
get_pairwise_corrs<-function(x,...){
  n = ncol(x);cn = colnames(x)
  m = diag(1,n)
  rownames(m) = cn; colnames(m) = cn
  for(i in 2:n){
    for(j in 1:(i-1)){
      inds = !apply(is.na(x[,c(i,j)]),1,any) & !apply(is.nan(x[,c(i,j)]),1,any)
      corr = cor(x[inds,i],x[inds,j],...)
      m[i,j] = corr; m[j,i]=corr
    }
  }
  return(m)
}

# Analyze the workload and time patterns
# Step 1: get the "hockey" sticks
# We take the exercise phase and plot the workloads vs time
# We want to know if we see big deviations from the expected 
# experiment in which the workloads manifest a roughly increasing
# function.
# ADD documentation (parameter names to their meaning)
get_segment_input <- function(row,val_mat,val,mat2,eps=1e-5){
  if(row%%100 == 0){print(row)}
  start_index = as.integer(get_first_index_that_fits_value(val,val_mat[row,]))
  end_index = as.integer(get_last_index_that_fits_value(val,val_mat[row,]))
  # If there is no start/end, we will return NA for our calculated values.
  if(is.na(start_index) || start_index == end_index){
    return(list(values=c(),start_index=-1,end_index=-1,is_mono=F))
  }
  #Move the end index back until it is not NA.
  while((is.na(mat2[row,end_index])||as.numeric(mat2[row,end_index])==0) && end_index > start_index){
    end_index = end_index - 1
  }
  #Move up the starting index until it is not NA.
  while((is.na(mat2[row,start_index])||as.numeric(mat2[row,start_index])==0) && start_index < end_index){
    start_index = start_index + 1
  }
  #print(start_index);print(end_index)
  # If there is no start/end, we will return values that indicate so.
  if(is.na(start_index) || start_index == end_index){
    return(list(values=c(),start_index=-1,end_index=-1,is_mono=F))
  }
  #Identify the index at which the first difference is identified.
  #We account for the fact that sometimes there is missing data points
  #in the middle of the recorded data points.
  start_point = end_index
  first_point_value = as.numeric(mat2[row,start_index])
  for(index in start_index:(end_index-1)){
    curr_point <- as.numeric(mat2[row,index])
    if(is.na(curr_point)){next}
    score_difference = curr_point - first_point_value
    if (score_difference > eps){
      start_point = index
      break
    }
  }
  # get some summary statistics
  v = as.numeric(mat2[row,start_point:end_index])
  #print(table(is.na(v)))
  if(length(v)<2){return(list(values=c(),start_index=-1,end_index=-1,is_mono=F))}
  is_monotone = all(v[1:(length(v)-1)]<=v[2:length(v)]+eps)
  out_list = list()
  out_list[["values"]] = v
  out_list[["start_index"]] = start_point
  out_list[["end_index"]] = end_index
  out_list[["is_mono"]] = is_monotone
  return(out_list)
}

# tests
#m1 = matrix(rep(1,100),nrow=5,ncol=20)
#for(j in 1:ncol(m1)){m1[,j]=j}
#m2 = matrix(rep("Exercise",100),nrow=5,ncol=20)
#get_segment_input(1,val_mat=m2, val="Exercise", mat2 = m1)
#m1[,1:2] = NA ; m1[,15:20] = NA
#get_segment_input(1,val_mat=m2, val="Exercise", mat2 = m1)
#m2[,2:16] = "A"
#get_segment_input(1,val_mat=m2, val="A", mat2 = m1)
#m2 = matrix(rep("Exercise",100),nrow=5,ncol=20)
#m2[,2:9] = "A"
#get_segment_input(1,val_mat=m2, val="A", mat2 = m1)

# source("https://bioconductor.org/biocLite.R")
# biocLite("impute",lib=getwd())
try({library("impute",lib.loc='/home/davidama/R/x86_64-pc-linux-gnu-library/3.2/')})
try({library('impute')})

# This function runs linear regression on a dependent variable y
# and a set of covariates
# OPTIONAL filter - additional_subcategories_to_exclude - which covariate categories to exclude
# feature2subcategory - a mapping of covariates into categories
get_lm_residuals<-function(y,covs,use_categorical=T,max_num_classes=5,
                           feature_is_numeric,additional_subcategories_to_exclude=NULL,feature2subcategory=NULL){
  inds = !is.na(y)
  if(!is.null(additional_subcategories_to_exclude)){
    col_inds = !is.element(colnames(covs),set=additional_subcategories_to_exclude)
    covs = covs[,col_inds]
  }
  y = y[inds]; covs = covs[inds,]
  x = covs[,feature_is_numeric[colnames(covs)]]
  x = as.matrix(x)
  mode(x) = 'numeric'
  x_imp = impute.knn(x)$data
  
  if(use_categorical){
    x2 = covs[,!feature_is_numeric[colnames(covs)]]
    new_x2 = c()
    for (j in 1:ncol(x2)){
      fx = x2[,j]
      fx_table = table(as.character(fx))
      if(sum(is.na(fx))/length(fx) >= 0.2){next}    
      if(length(fx_table)>max_num_classes || length(fx_table)<2){next}
      if(!is.factor(fx)){fx = as.factor(fx)}
      options(na.action='na.pass')
      fx_mat = model.matrix(~fx+0,data=data.frame(fx))
      colnames(fx_mat) = gsub(colnames(fx_mat),pattern="^fx",perl=T,replace="")
      colnames(fx_mat) = paste(colnames(x2)[j],colnames(fx_mat),sep='_')
      fx_mat = fx_mat[,apply(fx_mat,2,sd,na.rm=T)>0]
      #if(all(is.na(fx_mat))){break}
      if(length(fx_mat)==0){next}
      new_x2 = cbind(new_x2,fx_mat)
      #print(j)
    }
    #dim(new_x2)
    #per_Nas = apply(is.na(new_x2),2,sum)/nrow(new_x2)
    #sort(per_Nas)
    x2_imp = impute.knn(new_x2)$data
    x_imp = cbind(x_imp,x2_imp)
  }
  lm_obj = lm(y~x_imp,na.action=na.exclude)
  return(list(lm_obj=lm_obj,inds=inds))
}




################### Changes to some MTS functions ######################
VARXpred2<-function (m1, newxt = NULL, hstep = 1, orig = 0) 
{
  zt = as.matrix(m1$data)
  xt = as.matrix(m1$xt)
  p = m1$aror
  m = m1$m
  Ph0 = as.matrix(m1$Ph0)
  Phi = as.matrix(m1$Phi)
  Sig = as.matrix(m1$Sigma)
  beta = as.matrix(m1$beta)
  include.m = m1$include.mean
  nT = dim(zt)[1]
  k = dim(zt)[2]
  dx = dim(xt)[2]
  se = NULL
  if (length(Ph0) < 1) 
    Ph0 = matrix(rep(0, k), k, 1)
  if (hstep < 1) 
    hstep = 1
  if (orig < 1) 
    orig = nT
  if (length(newxt) > 0) {
    if (!is.matrix(newxt)) 
      newxt = as.matrix(newxt)
    h1 = dim(newxt)[1]
    hstep = min(h1, hstep)
    nzt = as.matrix(zt[1:orig, ])
    xt = rbind(xt[1:orig, , drop = FALSE], newxt)
    for (i in 1:hstep) {
      tmp = Ph0
      ti = orig + i
      for (i in 1:p) {
        idx = (i - 1) * k
        tmp = tmp + Phi[, (idx + 1):(idx + k)] %*% matrix(nzt[ti - 
                                                                i, ], k, 1)
      }
      if (m > -1) {
        for (j in 0:m) {
          jdx = j * dx
          tmp = tmp + beta[, (jdx + 1):(jdx + dx)] %*% 
            matrix(xt[ti - j, ], dx, 1)
        }
      }
      nzt = rbind(nzt, c(tmp))
    }
    mm = VARpsi(Phi, lag = hstep)
    Si = Sig
    se = matrix(sqrt(diag(Si)), 1, k)
    if (hstep > 1) {
      for (i in 2:hstep) {
        idx = (i - 1) * k
        wk = as.matrix(mm$psi[, (idx + 1):(idx + k)])
        Si = Si + wk %*% Sig %*% t(wk)
        se1 = sqrt(diag(Si))
        se = rbind(se, se1)
      }
    }
    #cat("Prediction at origin: ", orig, "\n")
    #cat("Point forecasts (starting with step 1): ", "\n")
    #print(round(nzt[(orig + 1):(orig + hstep), ], 5))
    #cat("Corresponding standard errors: ", "\n")
    #print(round(se[1:hstep, ], 5))
    return(cbind(
      round(nzt[(orig + 1):(orig + hstep), ], 5),
      round(se[1:hstep, ], 5)
    ))
  }
  else {
    cat("Need new data for input variables!", "\n")
  }
}
