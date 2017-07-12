# The following script loads the phenotypic data
# We start by running some QA tests in order to characterize the experiment that was
# taken for each individual
# On sherlock use: module load R

try({setwd("/Users/davidhsu/Documents/ukbb")})
try({setwd("/Users/david/Desktop/ukbb/")})
try({setwd('/scratch/PI/euan/projects/ukbb/da_dh/')})
source("auxiliary_functions.R")

### Constants for the analysis below ###
########################################
MERGE_SAME_WORKLOADS = T
MIN_REGR_SIZE = 3
# Paramaters by which we slice the data
REPEATS = c("0","1")
CATEGORIES = c("Category 1","Category 2")

# INPUT for the entire script:
# This code gets an object called subject_cleaned_pheno_data: a matrix of subjects (rows)
# vs the metadata features (columns)
# Load the preprocessed pheno data
load("biobank_collated_filtered_pheno_data.RData")
###################################

# Inspect the subject categories, are there discrepancies?
category_matrix = subject_cleaned_pheno_data[,get_regex_cols(colnames(subject_cleaned_pheno_data),"category")]
is_same_category<-function(x){
  if(any(is.na(x))){return(T)}
  return (x[1]==x[2])
}
table(category_matrix[,2])
category_discs = apply(category_matrix,1,is_same_category)
table(category_discs)
# look at the false rows
category_discs_rows = category_matrix[!category_discs,]

# Split by the experiment (time point to repeat number)
tp_to_repeat = sapply(colnames(subject_cleaned_pheno_data),function(x)strsplit(x,split = "\\.")[[1]][2])
table(tp_to_repeat)

###############################################
###############################################
########## Time series extraction #############
###############################################
###############################################

# each entry in the list below represents a slice of the metadata
# of time series. The name is the repeat;category.
# The entries are:
# pheno_data = the raw data matrix of the current slice
# regex matrices = the relevant time series columns partitioned by their regex 
# tps_vs_wkls = a list with the time points vs workloads vs heart rates (for each subject)
# regr_inputs = the preprocessed input for the regression analysis of the exercise phase
# subject2test_status
# subject2sex
preprocessed_data_slices = list()
for (REPEAT in REPEATS){
  for (CATEGORY in CATEGORIES){
    curr_slice_name = paste(REPEAT,CATEGORY,sep=';')
    
    # Slice the data, first by the repeat and then by the others
    # (otherwise we get duplications, e.g., for category)
    pheno_data = subject_cleaned_pheno_data[,tp_to_repeat==REPEAT]
    # now we have one category columns
    subj2category = pheno_data[,get_regex_cols(colnames(pheno_data),"category")]
    pheno_data = pheno_data[grepl(CATEGORY,subj2category),]
    
    # Create a list of matrices that correspond to the data gathered from the ECG Test.
    # The matrices created are: time, workload, phase, and heart rate.
    regex_matrices = extract_regex_matrices(pheno_data)  
    print(sapply(regex_matrices,dim))
    gc()
    
    # Process the data to get the HR, WD, and TP for each subject
    tps_vs_wkls = lapply(1:nrow(regex_matrices[["time"]]),get_segment_input,val_mat=regex_matrices[["phase"]],
                         val="Exercise", mat2 = regex_matrices[["workload"]],eps=1e-5)
    names(tps_vs_wkls)=rownames(regex_matrices$time)
    
    # Extract regression input and save the objects
    regr_inputs = list()
    for (j in 1:length(tps_vs_wkls)){
      if(j%%500 ==0){print(j)}
      curr_input =  tps_vs_wkls[[j]]
      curr_subj = names(tps_vs_wkls)[j]
      regr_inputs[[curr_subj]] = extract_regression_input_from_tp_and_wklds(curr_input)
    }
    
    # Other important parameters for the ECG analysis
    status_col = get_regex_cols(colnames(pheno_data),"status of test")
    subject2test_status = as.character(pheno_data[,status_col])
    names(subject2test_status) = rownames(pheno_data)
    sex_col = get_regex_cols(colnames(pheno_data),"Sex")
    subject2sex = NULL
    if(!is.null(sex_col)){
      subject2sex = as.character(pheno_data[,sex_col])
      names(subject2sex) = rownames(pheno_data)
    }
    
    # Store the results
    preprocessed_data_slices[[curr_slice_name]] = list()
    preprocessed_data_slices[[curr_slice_name]][["pheno_data"]] = pheno_data
    preprocessed_data_slices[[curr_slice_name]][["regex_matrices"]] = regex_matrices
    preprocessed_data_slices[[curr_slice_name]][["tps_vs_wkls"]] = tps_vs_wkls
    preprocessed_data_slices[[curr_slice_name]][["regr_inputs"]] = regr_inputs
    preprocessed_data_slices[[curr_slice_name]][["subject2test_status"]] = subject2test_status
    preprocessed_data_slices[[curr_slice_name]][["subject2sex"]] = subject2sex
  }
}
save(preprocessed_data_slices,file="fitness_analysis_preprocessed_data_slices.RData")

###############################################
###############################################
#################### End ######################
###############################################
###############################################

###############################################
###############################################
############### Tests and plots################
###############################################
###############################################

rand_ind = sample(1:nrow(regex_matrices[["time"]]))[1]
plot(regex_matrices[["time"]][rand_ind,],regex_matrices[["workload"]][rand_ind,])
plot(regex_matrices[["time"]][rand_ind,1:112],regex_matrices[["heartrate"]][rand_ind,])
#table(c(regex_matrices[["phase"]]))

# # TODO: check the MANUAL annotation, what does it mean?
# is_MANUAL_phase = regex_matrices[["phase"]] == "MANUAL"
# subjects_with_MANUAL = rowSums(is_MANUAL_phase,na.rm=T) > 0
# regex_matrices[["phase"]][which(subjects_with_MANUAL)[1],1:50]
# table(rowSums(is_MANUAL_phase,na.rm=T))
# table(colSums(is_MANUAL_phase,na.rm=T))
# ########################################################
# Analyze the phases
# We transform the phases into numbers separated by a 1000:
# makes the plotting easier.
numeric_phase_mat = regex_matrices[["phase"]]
numeric_phase_mat[numeric_phase_mat == "Pretest"] = "0"
numeric_phase_mat[numeric_phase_mat == "Exercise"] = "1"
# Rethink what to do with MANUAL after the TODO above is completed
numeric_phase_mat[numeric_phase_mat == "MANUAL"] = "1"
numeric_phase_mat[numeric_phase_mat=="Rest"] = "2"
mode(numeric_phase_mat) = "numeric"

# # This test: We wanted see that the last time point at each stage makes sense
# # compared to the explanation of UKBB.
# num_phases = 0:2
# for (phase in num_phases){
#   curr_m = numeric_phase_mat 
#   lst_inds = apply(curr_m,1,get_last_index_that_fits_value,val=phase)
#   subj2last_timepoint = c()
#   for(j in 1:nrow(curr_m)){
#     if(is.na(lst_inds[j])){
#       subj2last_timepoint[j] = NA
#       next
#     }
#     subj2last_timepoint[j] = regex_matrices[["time"]][j,lst_inds[j]]
#   }
#   hist(as.numeric(subj2last_timepoint))
#   table(as.numeric(subj2last_timepoint))
#   break
# }

# Visualize the experiment of a random subject
rand_ind = sample(1:nrow(regex_matrices[["time"]]))[1]
numeric_timepoints = regex_matrices[["time"]]
mode(numeric_timepoints) = "numeric"
new_time_order = numeric_phase_mat*1000 + numeric_timepoints
plot(new_time_order[rand_ind,],regex_matrices[["workload"]][rand_ind,],
     xlab = "Time (seconds within phases)", ylab="Workload")
#plot(regex_matrices[["heartrate"]][rand_ind,],regex_matrices[["workload"]][rand_ind,1:112])

# Some basic tests for the tps_vs_wkls object
# table(sapply(tps_vs_wkls,class))
# table(sapply(tps_vs_wkls,function(x)class(x)=="list"&&x$is_mono)) / length(tps_vs_wkls)
# table(sapply(tps_vs_wkls,function(x)class(x)=="list"&&x$is_mono))
# num_wlds = sapply(tps_vs_wkls,function(x)length(x$values))
# not_mono = sapply(tps_vs_wkls,function(x)class(x)=="list"&&!x$is_mono)
# not_mono_with_wlds = not_mono & num_wlds>0
# names(subject2test_status) = names(not_mono)
# table(sapply(tps_vs_wkls[not_mono],function(x)length(x$values)))
# table(not_mono_with_wlds,subject2test_status[names(not_mono_with_wlds)])
# # look at the extracted window size
# hist(num_wlds)

###############################################
###############################################
#################### End ######################
###############################################
###############################################


##########################################################################################
load("regression_analysis_input_data.RData")
# Changing this to regr_inputs means running 
# the regression analyses below on the unmerged data
lm_inputs = regr_inputs_wd_merged
subject2test_status = subject2test_status[names(lm_inputs)]
subject2sex = subject2sex[names(lm_inputs)]

# Our main goal below is to get an estimation for the subjects' fitness.
# Preprocessing: we use the subject status, window size, mean/max HR and mean/max WD
sizes = sapply(lm_inputs,function(x)length(x$HR))
# The warnings we get from the commands below are
# due to subjects with no data
mean_HRs = sapply(lm_inputs,function(x)mean(x$HR))
max_HRs = sapply(lm_inputs,function(x)max(x$HR))
sd_HRs = sapply(lm_inputs,function(x)sd(x$HR))
mean_WDs = sapply(lm_inputs,function(x)mean(x$WD))
max_WDs = sapply(lm_inputs,function(x)max(x$WD))
has_completed = subject2test_status == "Fully completed"
has_completed = has_completed[names(sizes)]
table(has_completed)
has_completed_and_has_small_window = has_completed & sizes < 10
#table(has_completed_and_has_small_window)
subject_technical_class = rep("C:Large",length(sizes))
names(subject_technical_class) = names(sizes)
subject_technical_class[!has_completed] = "Not C"
subject_technical_class[has_completed_and_has_small_window] = "C:Small"
subject_technical_class[sizes<1] = "No data"
table(subject_technical_class)
fully_comps = names(which(subject_technical_class=="C:Large"))
max_WD_pheno_data = as.numeric(pheno_data[,"Maximum workload during fitness test.0.0"])
names(max_WD_pheno_data) = rownames(pheno_data)
all(names(max_WD_pheno_data)==names(max_WDs))
# all the infinite cases are when our WD is null since no increasing
# window was found. In these cases, we take the max WD directly from the
# pheno data
# # QA
# inds = which(sapply(lm_inputs,function(x)any(is.infinite(max(x$WD)))))
# all(sapply(lm_inputs[inds],function(x)is.null(x$WD)))
max_WDs[is.infinite((max_WDs))] =  max_WD_pheno_data[is.infinite((max_WDs))]
plot(max_WD_pheno_data,max_WDs)

# # QA: check that all names are the same as the regex matrices
# length(sizes)==nrow(regex_matrices$time)
# all(names(sizes)==rownames(regex_matrices$time))
# length(subject_technical_class)==nrow(regex_matrices$time)
# all(names(subject_technical_class)==rownames(regex_matrices$time))
# length(mean_HRs)==nrow(regex_matrices$time)
# all(names(mean_HRs)==rownames(regex_matrices$time))
# length(mean_WDs)==nrow(regex_matrices$time)
# all(names(mean_WDs)==rownames(regex_matrices$time))

# Get some boxplots to show statistics based on completion status
# We differentiate between subjects that completed the exercise stage
# and had a relatively large window size (of increasing WLDs) and those
# with a small one (<10). 
par(mfrow=c(2,3),mar=c(4,4,4,4))
# Warning in the commands below are due to subjects with no data, ignore
boxplot(get_list_by_values(sizes,subject_technical_class),las=2,ylab="window size")
boxplot(get_list_by_values(mean_WDs,subject_technical_class),las=2,ylab = "mean workload")
boxplot(get_list_by_values(max_WDs,subject_technical_class),las=2,ylab="max workload")
boxplot(get_list_by_values(mean_HRs,subject_technical_class),las=2,ylab="mean heart rate")
boxplot(get_list_by_values(max_HRs,subject_technical_class),las=2,ylab="max heart rate")
boxplot(get_list_by_values(sd_HRs,subject_technical_class),las=2,ylab="heart rate sd")
# Conclusions:
# Subjects that did not complete the tests (most likely have) low fitness.
# Subject that completed the test and had very few WDs are more similar to the 
# uncompleted cases than to the other completed cases.
# Changing the window size threshold for 15 eliminates the separation between the completed
# cases in terms of hear rates (mean, max).

# # Subject classes vs. sex
# tbl = table(subject_technical_class,subject2sex[names(subject_technical_class)])
# tbl/colSums(tbl)
# chisq.test(table(subject_technical_class,subject2sex[names(subject_technical_class)]))$p.value

# Sanity checks: correlation with time of the exercise and rest phases
rho_exercise = c();wsize_exercise = c()
rho_rest = c();wsize_rest=c()
M_ph = regex_matrices$phase[,1:112]
M_T = regex_matrices$time[,1:112]
M_hr = regex_matrices$heartrate
mode(M_hr) = 'numeric'
mode(M_T) = 'numeric'
corMethod="spearman"
for(nn in rownames(regex_matrices$time)){
  inds1 = M_ph[nn,]=="Exercise"
  inds1 = inds1 & !is.na(M_hr[nn,])
  rho_exercise[nn] = cor(M_T[nn,inds1],M_hr[nn,inds1],method=corMethod)
  wsize_exercise[nn] = sum(inds1)
  inds2 = M_ph[nn,]=="Rest"
  inds2 = inds2 & !is.na(M_hr[nn,])
  rho_rest[nn] = cor(M_T[nn,inds2],M_hr[nn,inds2],method=corMethod)
  wsize_rest[nn] = sum(inds2)
}
table(rho_exercise>0.5)/sum(!is.na(rho_exercise))
table(rho_rest< -0.5)/sum(!is.na(rho_rest))
# plot(pairwise_plot(rho_exercise,rho_rest,hexbin,fully_comps,xbins=100),xlab="Rho exercise",ylab="Rho rest")
# plot(pairwise_plot(wsize_exercise,rho_exercise,hexbin,fully_comps,xbins=100),xlab="Window size",ylab="Rho",main="Exercise")
# boxplot(get_list_by_values(rho_rest,wsize_rest,5),xlab="Window size",ylab="Rho",main="Exercise")
# par(mfrow=c(2,2))
# l = get_list_by_values(rho_exercise,subject_technical_class,values_to_ignore = "No data")
# boxplot(l,ylab="Rho",main="Exercise",las=2)
# qqplot(l[[1]],l[[3]],xlab=names(l)[1],ylab=names(l)[3],cex=1.3,col="grey",pch=20,main="Exercise rho");abline(0,1,lty=2,lwd=2)
# l = get_list_by_values(rho_rest,subject_technical_class,values_to_ignore = "No data")
# boxplot(l,ylab="Rho",main="Rest",las=2)
# qqplot(l[[1]],l[[3]],xlab=names(l)[1],ylab=names(l)[3],cex=1.3,col="grey",pch=20,main="Rest rho");abline(0,1,lty=2,lwd=2)

# Compare the R-sq scores of different regression analyses
r2_values = c()
for(i in names(lm_inputs)){
  if(is.null(lm_inputs[[i]])){next}
  WD = lm_inputs[[i]]$WD; HR = lm_inputs[[i]]$HR
  TP = lm_inputs[[i]]$TP
  if(length(WD)<3){
    r2_values = rbind(r2_values,c(NA,NA,NA))
    rownames(r2_values)[nrow(r2_values)] = i
    next
  }
  obj1 = lm(HR~WD); obj2 = lm(HR~TP); obj3 = lm(HR~TP+WD)
  r2_values = rbind(r2_values,c(get_lm_r2(obj1),get_lm_r2(obj2),get_lm_r2(obj2)))
  rownames(r2_values)[nrow(r2_values)] = i
  if(as.numeric(i)%%1000==0){print(i)}
}

# Compute  the autocorrelations of the heart rates
autocorrs = c()
for(i in names(lm_inputs)){
  if(is.null(lm_inputs[[i]])){next}
  WD = lm_inputs[[i]]$WD; HR = lm_inputs[[i]]$HR
  TP = lm_inputs[[i]]$TP
  if(length(WD) < 5){
    autocorrs = rbind(autocorrs,c(NA,NA,NA))
    rownames(autocorrs)[nrow(autocorrs)] = i
    next
  }
  corrs = acf(HR,plot=F,lag.max=3,type="correlation")$acf
  autocorrs = rbind(autocorrs,as.numeric(corrs)[-1])
  rownames(autocorrs)[nrow(autocorrs)] = i
  if(as.numeric(i)%%1000==0){print(i)}
}


try({
  r2 = r2_values[,1]
  autoc = autocorrs[,1]
  library(hexbin)
  # pariwise plots
  plot(pairwise_plot(r2,autoc,hexbin,fully_comps,xbins=100),ylab="Autocorrelation",xlab="R^2")
  #plot(pairwise_plot(r2,wd_tp_cor,hexbin,fully_comps,xbins=100),ylab="Dev",xlab="R^2")
  # boxplots
  par(mfrow=c(1,2),mar=c(5,4,4,4))
  boxplot(get_list_by_values(r2,subject_technical_class[names(r2)]),ylim=c(0.6,1),ylab="R^2",las=2)
  boxplot(get_list_by_values(autocorrs[,1],subject_technical_class[rownames(autocorrs)]),las=2,ylab="Autocorrelation")
  plot(pairwise_plot(r2,autoc,hexbin,fully_comps,xbins=100),ylab="Autocorrelation",xlab="R^2")
  plot(pairwise_plot(r2,sd_HRs,hexbin,fully_comps,xbins=100),ylab="HR sd",xlab="R^2")
})

# Save the scores
save(rho_exercise,rho_rest,r2_values,autocorrs,file="Regression_correlation_time_series_qc_scores.RData")


#################################### Compute the fitness scores #####################################
load("Regression_correlation_time_series_qc_scores.RData")
load("UKBB_phenotypic_data_for_GWAS.RData")
# Get predictions for specific WDs during the exercise phase
# NA: subjects with low R2
# -1: subjects not in the C:Large class
subject_ols_preds_50 = c();subject_ols_preds_100 = c();subject_ols_preds_150 = c()
lm_objs = list()
for(i in names(lm_inputs)){
  # if (subject_technical_class[i]!="C:Large"){
  #   subject_ols_preds_50[i]=-1
  #   subject_ols_preds_100[i]=-1
  #   subject_ols_preds_150[i]=-1
  #   next
  # }
  if(is.null(lm_inputs[[i]]) || is.na(rho_exercise[i]) || rho_exercise[i] < 0.5){
    subject_ols_preds_50[i]=NA
    subject_ols_preds_100[i]=NA
    subject_ols_preds_150[i]=NA
    next
  }
  WD = lm_inputs[[i]]$WD; HR = lm_inputs[[i]]$HR;TP = lm_inputs[[i]]$TP
  if(is.null(WD) || length(WD) < 5){
    subject_ols_preds_50[i]=NA
    subject_ols_preds_100[i]=NA
    subject_ols_preds_150[i]=NA
    next
  }
  lm_objs[[i]] = lm(HR~WD)
  if (lm_objs[[i]]$coefficients[2]<0 | get_lm_r2(lm_objs[[i]]) < 0.5){
    subject_ols_preds_50[i]=NA
    subject_ols_preds_100[i]=NA
    subject_ols_preds_150[i]=NA
    next
  }
  subject_ols_preds_50[i] = predict(lm_objs[[i]],data.frame(WD=50))
  subject_ols_preds_100[i] = predict(lm_objs[[i]],data.frame(WD=100))
  subject_ols_preds_150[i] = predict(lm_objs[[i]],data.frame(WD=150))
  if(as.numeric(i)%%1000==0){print(i)}
}
all(names(subject_ols_preds_50)==rownames(regex_matrices[[1]]))
inds = subject_ols_preds_150!=-1 & !is.na(subject_ols_preds_150)
cor(subject_ols_preds_150[inds],subject_ols_preds_50[inds])

# # Look at some strong negative correlations
# table(sapply(lm_objs,function(x)x$coefficients[2])<0)
# neg_corrs = names(which(sapply(lm_objs,function(x)unname(x$coefficients[2]))<0))
# i=neg_corrs[1]
# WD = lm_inputs[[i]]$WD; HR = lm_inputs[[i]]$HR;TP = lm_inputs[[i]]$TP
# plot(WD,HR)
#############################################################

# Get an estimation of the decline in HR, relative to WD during the resting phase
subj2rest_data = list()
for (nn in rownames(regex_matrices$phase)){
  inds = regex_matrices$phase[nn,] == "Rest"
  inds = inds[1:112]
  inds[is.na(regex_matrices$heartrate[nn,]) | is.na(inds)] = F
  rest_WD = regex_matrices$workload[nn,inds]
  rest_HR = regex_matrices$heartrate[nn,inds]
  rest_TP = as.numeric(regex_matrices$time[nn,inds])
  subj2rest_data[[nn]] = list(WD=rest_WD,HR=rest_HR,TP=rest_TP)
  if(as.numeric(nn)%%1000==0){print(nn)}
}
#sort(table(unlist(sapply(subj2rest_data,function(x)x$TP))))

rest_HR_decline_score = c()
rest_HR_decline_WD_normalized_score = c()
rest_HR_decline_HR_normalized_score = c()
max_time_diff = 3
rest_times = c(29,60)
for(nn in names(subject_ols_preds_150)){
  TP = subj2rest_data[[nn]]$TP
  HR = as.numeric(subj2rest_data[[nn]]$HR)
  WD = as.numeric(subj2rest_data[[nn]]$WD)
  v1=c();v2=c();v3=c()
  if(is.na(rho_rest[nn]) || rho_rest[nn]> -0.5){
    rest_HR_decline_score = rbind(rest_HR_decline_score,rep(NA,length(rest_times)))
    rest_HR_decline_WD_normalized_score = rbind(rest_HR_decline_WD_normalized_score,rep(NA,length(rest_times)))
    rest_HR_decline_HR_normalized_score = rbind(rest_HR_decline_HR_normalized_score,rep(NA,length(rest_times)))
    next
  }
  try({
    HR_smoothed = smooth.spline(HR)
    #par(mfrow=c(1,2));plot(HR);plot(HR_smoothed)
    HR = HR_smoothed$y
  })
  for(tt in rest_times){
    tt_char = as.character(tt)
    t_diffs = abs(TP-tt)
    inds = which(t_diffs<=max_time_diff)
    if(length(inds)==0 || sum(inds)==0){
      v1[tt_char]=NA;v2[tt_char]=NA;v3[tt_char]=NA;next
    }
    tt_ind = which(t_diffs==min(t_diffs))[1]
    v1[tt_char] = HR[1] - HR[tt_ind]
    v2[tt_char] = v1[tt_char]*max_WDs[nn]
    v3[tt_char] = HR[tt_ind]/HR[1]
  }
  rest_HR_decline_score = rbind(rest_HR_decline_score,v1)
  rest_HR_decline_WD_normalized_score = rbind(rest_HR_decline_WD_normalized_score,v2)
  rest_HR_decline_HR_normalized_score = rbind(rest_HR_decline_HR_normalized_score,v3)
  if(as.numeric(nn)%%1000==0){print(nn)}
}
rownames(rest_HR_decline_WD_normalized_score) = names(subject_ols_preds_150)
rownames(rest_HR_decline_HR_normalized_score) = names(subject_ols_preds_150)
rownames(rest_HR_decline_score) = rownames(rest_HR_decline_WD_normalized_score)
colnames(rest_HR_decline_WD_normalized_score) = paste("HR_rest_diff_wd",rest_times)
colnames(rest_HR_decline_HR_normalized_score) = paste("HR_rest_diff_perc",rest_times)
colnames(rest_HR_decline_score) = paste("HR_rest_diff",rest_times)

# Some QA
#all(names(max_WDs)==names(subject_ols_preds_50))
#length(names(max_WDs)) == nrow(regex_matrices$time)

# Combine all fitness scores
fitness_scores = cbind(
                       subject_ols_preds_50,
                       subject_ols_preds_100,
                       subject_ols_preds_150,
                       rest_HR_decline_score,
                       rest_HR_decline_WD_normalized_score,
                       rest_HR_decline_HR_normalized_score,
                       max_WD_pheno_data
                  )
colnames(fitness_scores)[1:3] = gsub(colnames(fitness_scores)[1:3],pattern="subject_",replace="")
colnames(fitness_scores)[1:3] = gsub(colnames(fitness_scores)[1:3],pattern="ols_preds",replace="HR_pred")
save(fitness_scores,pheno_data,file="UKBB_phenotypic_data_for_GWAS.RData")
