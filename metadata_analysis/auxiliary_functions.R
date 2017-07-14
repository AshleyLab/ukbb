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
get_lm_residuals<-function(y,covs,use_categorical=T,max_num_classes=5, max_allowed_na_per = 0.2,
                           feature_is_numeric,additional_subcategories_to_exclude=NULL){
  inds = !is.na(y)
  if(!is.null(additional_subcategories_to_exclude)){
    col_inds = !is.element(colnames(covs),set=additional_subcategories_to_exclude)
    covs = covs[,col_inds]
  }
  y = y[inds]
  covs = covs[names(y),]
  x = covs[,feature_is_numeric[colnames(covs)]]
  x = as.matrix(x)
  mode(x) = 'numeric'
  print ("Imputing missing values in numeric part")
  x_imp = impute.knn(x)$data
  print ("Done imputing missing values in numeric part")
  
  if(use_categorical && sum(!feature_is_numeric[colnames(covs)])>0){
    x2 = covs[,!feature_is_numeric[colnames(covs)]]
    new_x2 = c()
    for (j in 1:ncol(x2)){
      fx = x2[,j]
      fx_table = table(as.character(fx))
      if(sum(is.na(fx))/length(fx) >= max_allowed_na_per){next}    
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
    }
    x2_imp = impute.knn(new_x2)$data
    x_imp = cbind(x_imp,x2_imp)
  }
  print("Done creating the covariate matrix for the regression analysis")
  print(dim(x_imp))
  lm_obj = lm(y~x_imp,na.action=na.exclude)
  return(list(lm_obj=lm_obj,inds=inds))
}

# Create a list of matrices that correspond to the data gathered from the ECG Test.
# The matrices created are: time, workload, phase, and heart rate.
# This function gets a pheno data matrix that contains a single time series
# (at most) per subject.
extract_regex_matrices<-function(pheno_data){
  regex_list = list()
  regex_list[["time"]] = "ECG, phase time"
  regex_list[["workload"]] = "ECG, load"
  regex_list[["phase"]] = "ECG, trend phase name"
  regex_list[["heartrate"]] = "ECG, heart rate"
  regex_cols = list()
  for(nn in names(regex_list)){
    regex_cols[[nn]] = get_regex_cols(colnames(pheno_data),regex_list[[nn]])
  }
  regex_matrices = list()
  for(nn in names(regex_list)){
    regex_matrices[[nn]] = as.matrix(pheno_data[,regex_cols[[nn]]])
  }
  return(regex_matrices)
}

# The input here is an entry of tps_vs_wkls.
# These objects are obtained after running get_segment_input
extract_regression_input_from_tp_and_wklds<-function(curr_input){
  if(!curr_input$is_mono || length(curr_input$value)<MIN_REGR_SIZE){
    return(list(HR=c(),WD=c(),TP=c()))
  }
  ind1 = curr_input$start_index
  ind2 = curr_input$end_index
  x = as.numeric(regex_matrices[["heartrate"]][j,ind1:ind2])
  y = as.numeric(regex_matrices[["workload"]][j,ind1:ind2])
  z = as.numeric(regex_matrices$time[j,ind1:ind2])
  inds = !is.na(x) & !is.na(y)
  x=x[inds];y=y[inds];z=z[inds]
  if(length(z)<MIN_REGR_SIZE){
    return(list(HR=c(),WD=c(),TP=c()))
  }
  if(MERGE_SAME_WORKLOADS){
    HR = as.numeric(tapply(x,y,mean,na.rm=T))
    TP = as.numeric(tapply(z,y,mean,na.rm=T))
    WD = as.numeric(tapply(y,y,mean,na.rm=T))
    return(list(HR=HR,WD=WD,TP=TP))
  }
  else{
    return(list(HR=x,WD=y,TP=z))
  }
}

# This function first performs a set of tests:
# 1. That the quality rho_exercise score is high
# 2. That there are time points (i.e., > MIN_REGR_SIZE)
# 3. That the regression coefficient is positive and the R2 is high
# If the current subject data fails in any of them
# then the functions returns NA
# Otherwise an lm object is returned.
get_lm_object_for_hr_prediction<-function(lm_input,rho_exercise,
                                          thr1=RHO_Q_THRESHOLD,thr2=R2_Q_THRESHOLD){
  if(is.null(lm_input) || is.na(rho_exercise) || rho_exercise< thr1){
    return(NA)
  }
  WD = lm_input$WD; HR = lm_input$HR;TP = lm_input$TP
  if(is.null(WD) || length(WD) < 5){
    return(NA)
  }
  lm_obj = lm(HR~WD)
  if (lm_obj$coefficients[2]<0 | get_lm_r2(lm_obj) < thr2){
    return(NA)
  }
  return (lm_obj)
}

get_prediction_from_lm_object_for_hr<-function(obj,WD){
  if(is.null(obj) || is.na(obj)){
    return(NA)
  }
  return (unname(predict(obj,data.frame(WD=WD))))
}

# Like the exercise phase analysis, this function tests
# the data and then computes the score.
# If any of the following tests fail, the function returns NA:
# 1. There is a rho_test score and it is high
# 2. There is a time point that is close enough to the requested one
get_rest_phase_fitness_score<-function(rest_data,rho_rest,rest_time = 60,
                                       smooth_data=T,use_ratio=T,max_time_diff = 3){
  TP = rest_data$TP
  HR = as.numeric(rest_data$HR)
  WD = as.numeric(rest_data$WD)
  v1=c();v2=c();v3=c()
  if(is.na(rho_rest) || rho_rest> -0.5){
    return(NA)
  }
  if(smooth_data && length(HR)>3){
    try({
      HR_smoothed = smooth.spline(HR)
      HR = HR_smoothed$y
    })
  }
  tt = rest_time
  tt_char = as.character(tt)
  t_diffs = abs(TP-tt)
  inds = which(t_diffs<=max_time_diff)
  if(length(inds)==0 || sum(inds)==0){return(NA)}
  tt_ind = which(t_diffs==min(t_diffs))[1]
  if(!use_ratio){return(HR[1] - HR[tt_ind])}
  return(HR[tt_ind]/HR[1])
}

# A function that merges the results of the fitness scores computation
merge_scores_list<-function(l,plot_class_averages = T,...){
  all_subjects = c()
  for(x in l){
    all_subjects = union(all_subjects,names(x))
  }
  l = lapply(l,function(x)x[!is.na(x)])
  scores = rep(0,length(all_subjects));classes = rep(NA,length(all_subjects));counts = rep(0,length(all_subjects))
  names(scores) = all_subjects;names(classes) = all_subjects;names(counts) = all_subjects
  for(nn in names(l)){
    curr_scores = l[[nn]]
    curr_names = names(curr_scores)
    scores[curr_names] = scores[curr_names] + curr_scores
    counts[curr_names] = counts[curr_names] + 1
    curr_classes = classes[curr_names]
    curr_na_classes = is.na(curr_classes)
    curr_classes[curr_na_classes] = nn
    curr_classes[!curr_na_classes] = "multiple"
    classes[curr_names] = curr_classes
  }
  table(counts)
  tests = all(is.na(classes[counts==0]))
  tests = tests & all(names(which(counts=="multiple")) == names(which(counts>1)))
  scores[counts>1] = scores[counts>1]/counts[counts>1]
  scores[counts==0] = NA
  d = data.frame(scores=scores,classes=classes)
  if(plot_class_averages){plot.design(d,...)}
  return(d)
}

