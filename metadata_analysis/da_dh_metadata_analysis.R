# The following script loads the phenotypic data
# We start by running some QA tests in order to characterize the experiment that was
# taken for each individual

# On sherlock use: module load R

try({setwd("/Users/davidhsu/Documents/ukbb")})
try({setwd("/Users/david/Desktop/ukbb/")})
try({setwd('/scratch/PI/euan/projects/ukbb/da_dh/')})
# We cannot change the working directory as we work in RStudio, so we comment out this line temporarily.
source("auxiliary_functions.R")

# ####################################################################################################################
# ####################################################################################################################
# Load the following data structures and functions to start the workspace
# Load the preprocessed pheno data
load("biobank_collated_filtered_pheno_data.RData")

# Split by the experiment
tp_to_repeat = sapply(colnames(subject_cleaned_pheno_data),function(x)strsplit(x,split = "\\.")[[1]][2])

# Paramaters by which we slice the data
REPEAT = "0"
CATEGORY = "Category 1"

# Slice the data, first by the repeat and then by the others
# (otherwise we get duplications, e.g., for category)
pheno_data = subject_cleaned_pheno_data[,tp_to_repeat==REPEAT]
# now we have one category columns
subj2category = pheno_data[,get_regex_cols(colnames(pheno_data),"category")]
pheno_data = pheno_data[grepl(CATEGORY,subj2category),]
# remove uneeded objects for the subsequent analyses
rm(subject_cleaned_pheno_data);rm(tp_to_repeat);rm(subj2category)
dim(pheno_data)
gc()

#Create a list of matrices that correspond to the data gathered from the ECG Test.
# The matrices created are time, workload, phase, and heart rate.
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
sapply(regex_matrices,dim)
gc()

# Other very important parameters for the ECG analysis
status_col = get_regex_cols(colnames(pheno_data),"status of test")
subject2test_status = as.character(pheno_data[,status_col])
names(subject2test_status) = rownames(pheno_data)
table(subject2test_status)
sex_col = get_regex_cols(colnames(pheno_data),"Sex")
subject2sex = as.character(pheno_data[,sex_col])
names(subject2sex) = rownames(pheno_data)
table(subject2sex)

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

##################### Analyze the phases #####################
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

tps_vs_wkls = lapply(1:nrow(regex_matrices[["time"]]),get_segment_input,val_mat=regex_matrices[["phase"]],
                         val="Exercise", mat2 = regex_matrices[["workload"]],eps=1e-5)
names(tps_vs_wkls)=rownames(regex_matrices$time)

# Some basic tests
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

# Extract regression input and save the objects
MIN_REGR_SIZE = 3
regr_inputs = list()
regr_inputs_wd_merged = list()
for (j in 1:length(tps_vs_wkls)){
  print(j)
  curr_input =  tps_vs_wkls[[j]]
  curr_subj = names(tps_vs_wkls)[j]
  if(!curr_input$is_mono || length(curr_input$value)<MIN_REGR_SIZE){
    regr_inputs[[curr_subj]] = list(HR=c(),WD=c(),TP=c())
    regr_inputs_wd_merged[[curr_subj]] = list(HR=c(),WD=c(),TP=c())
    next
  }
  ind1 = curr_input$start_index
  ind2 = curr_input$end_index
  x = as.numeric(regex_matrices[["heartrate"]][j,ind1:ind2])
  y = as.numeric(regex_matrices[["workload"]][j,ind1:ind2])
  z = as.numeric(regex_matrices$time[j,ind1:ind2])
  inds = !is.na(x) & !is.na(y)
  x=x[inds];y=y[inds];z=z[inds]
  if(length(z)<MIN_REGR_SIZE){
    regr_inputs[[curr_subj]] = list(HR=c(),WD=c(),TP=c())
    regr_inputs_wd_merged[[curr_subj]] = list(HR=c(),WD=c(),TP=c())
    next
  }
  regr_inputs[[curr_subj]] = list(HR=x,WD=y,TP=z)
  HR = as.numeric(tapply(x,y,mean,na.rm=T))
  TP = as.numeric(tapply(z,y,mean,na.rm=T))
  WD = as.numeric(tapply(y,y,mean,na.rm=T))
  regr_inputs_wd_merged[[curr_subj]] = list(HR=HR,WD=WD,TP=TP)
}
save(tps_vs_wkls,regr_inputs,regr_inputs_wd_merged,regex_matrices,
     subject2test_status,subject2sex,
     file="regression_analysis_input_data.RData")

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

####### Get a bunch of covariates and regress vis a vis our scores ###############
load("UKBB_phenotypic_data_for_GWAS.RData")
# Check which features are numeric
feature_is_numeric = c()
for(j in 1:ncol(pheno_data)){
  feature_is_numeric[j] = !all(is.na(as.numeric(as.character(pheno_data[,j]))))
}
#table(feature_is_numeric)
names(feature_is_numeric) = colnames(pheno_data)
bin_features = apply(pheno_data,2,function(x)length(unique(x))==2)
#names(which(bin_features))
# Remove features with too many NAs and those that are ECG-related
# Exclude features that are 20% or more NAs
all_features = colnames(pheno_data)
features_to_exclude = apply(is.na(pheno_data),2,sum)/nrow(pheno_data) >= 0.2
for(j in 1:ncol(pheno_data)){
  if (features_to_exclude[j]){next}
  fx = as.numeric(as.character(pheno_data[,j]))
  if(feature_is_numeric[j] && sum(is.na(fx))/length(fx) >= 0.2){
    features_to_exclude[j]=T
  }
}
selected_features = all_features[!features_to_exclude]
selected_features = selected_features[!grepl(selected_features,pattern="ECG")]
initial_cov_matrix = pheno_data[,selected_features]
initial_cov_matrix = initial_cov_matrix[,apply(initial_cov_matrix,2,function(x)length(unique(x[!is.na(x)]))) >1]
dim(initial_cov_matrix)
table(feature_is_numeric[colnames(initial_cov_matrix)])
#######################
# Load mapping of feature names to categories
feature_metadata =  read.csv('Data_Dictionary_Showcase.csv')
library(xlsx)
feature_data_recommentadions = read.xlsx2("HR_exercise_pred_100_covariate_analysis.xlsx",1)
feature_recommentadions = as.character(feature_data_recommentadions[,ncol(feature_data_recommentadions)])
names(feature_recommentadions) = feature_data_recommentadions[,1]
feature_recommentadions = feature_recommentadions[names(feature_recommentadions)!=""]

#feature_metadata =  read.csv('../Data_Dictionary_Showcase.csv')
feature_category_data = as.character(feature_metadata[,"Cat2_Title"])
names(feature_category_data) = as.character(feature_metadata[,"Field"])
feature_subcategory_data = as.character(feature_metadata[,"Cat3_Title"])
names(feature_subcategory_data) = as.character(feature_metadata[,"Field"])
initial_cov_matrix_feature_category = feature_category_data[sapply(colnames(initial_cov_matrix),function(x)strsplit(x,split='\\.')[[1]][1])]
names(initial_cov_matrix_feature_category) = colnames(initial_cov_matrix)
initial_cov_matrix_feature_subcategory = feature_subcategory_data[sapply(colnames(initial_cov_matrix),function(x)strsplit(x,split='\\.')[[1]][1])]
names(initial_cov_matrix_feature_subcategory) = colnames(initial_cov_matrix)
# table(initial_cov_matrix_feature_category)
# table(initial_cov_matrix_feature_subcategory)
# table(initial_cov_matrix_feature_category,initial_cov_matrix_feature_subcategory)
# table(is.na(initial_cov_matrix_feature_category))
initial_cov_matrix_is_numeric = feature_is_numeric[colnames(initial_cov_matrix)]
table(initial_cov_matrix_feature_category,initial_cov_matrix_is_numeric)
# Some manual inspection:
initial_cov_matrix_feature_category[get_regex_cols(colnames(initial_cov_matrix),"blood",ignore.case=T)]
initial_cov_matrix_feature_category["Target heart rate achieved.0.0"]
# Some of these are relevant and some are not -  a mixture
initial_cov_matrix_feature_category[initial_cov_matrix_feature_category == "Physical measures"]
# Blood measures is in here
initial_cov_matrix_feature_category[initial_cov_matrix_feature_category == "Assay results"]
# Others
initial_cov_matrix_feature_subcategory[initial_cov_matrix_feature_subcategory == "Baseline characteristics"]
feature_is_numeric[names(which(initial_cov_matrix_feature_subcategory == "Baseline characteristics"))]
initial_cov_matrix_feature_category[initial_cov_matrix_feature_category == "Cognitive function"]
table(initial_cov_matrix_feature_category)
table(initial_cov_matrix_feature_subcategory)
initial_cov_matrix_feature_subcategory[initial_cov_matrix_feature_subcategory=="DXA assessment"]
initial_cov_matrix_feature_subcategory[initial_cov_matrix_feature_subcategory=="Anthropometry"]
hist(initial_cov_matrix[,"Basal metabolic rate.0.0"])
plot(fitness_scores[,2],initial_cov_matrix[,"Basal metabolic rate.0.0"])

# Look at medication and job
cols = colnames(initial_cov_matrix)
fs = cols[grepl("medic",cols,ignore.case = T)]
apply(initial_cov_matrix[,fs],2,table)

categories_to_exclude = c("Cognitive function")
potential_subcategories_to_exclude = c()
subcategories_to_exclude = c("ECG during exercise")
features_to_exclude = names(feature_recommentadions)[grepl("ignore",feature_recommentadions,ignore.case = T)
                                                     | grepl("score",feature_recommentadions,ignore.case = T)]
for(j in 1:ncol(initial_cov_matrix)){
  currf = colnames(initial_cov_matrix)[j]
  if(is.element(initial_cov_matrix_feature_category[currf],set=categories_to_exclude)){
    features_to_exclude = c(features_to_exclude,currf)
    next
  }
  if(is.element(initial_cov_matrix_feature_subcategory[currf],set=subcategories_to_exclude)){
    features_to_exclude = c(features_to_exclude,currf)
    next
  }
  # Final test: if feature is non-numeric and breakes the samples into
  # classes that are too small
  if(initial_cov_matrix_is_numeric[j]){next}
  fcounts = table(initial_cov_matrix[,j])
  if(max(fcounts)<500){
    features_to_exclude = c(features_to_exclude,colnames(initial_cov_matrix)[j])
    print("excluding")
  }
}
length(features_to_exclude)
covariate_matrix = initial_cov_matrix[,!is.element(colnames(initial_cov_matrix),set=features_to_exclude)]
dim(covariate_matrix)
# table(initial_cov_matrix_feature_category[colnames(covariate_matrix)])
# table(initial_cov_matrix_feature_subcategory[colnames(covariate_matrix)])

# # Look at co-linearity among some features and run PCA
# pca_matrices = list()
# m = pheno_data_cov_submatrices[["fat_features"]]
# rnames = rownames(pheno_data_cov_submatrices[["fat_features"]])
# m = as.matrix(m)
# m = m[(apply(is.na(m),1,sum))<5,]
# m = impute.knn(m)$data
# pca_m = prcomp(m,scale. = T)
# plot(pca_m)
# PCs = pca_m$x[,1:2]
# pca_m$rotation[,1:2]
# corrplot(cor(m))
# missing_instances = matrix(NA,nrow = length(setdiff(rnames,rownames(PCs))),ncol=2)
# rownames(missing_instances) = setdiff(rnames,rownames(PCs))
# PCs = rbind(PCs,missing_instances)
# PCs = PCs[rnames,]
# colnames(PCs) = paste("fat",colnames(PCs))

fitness_scores_inds = c(2,9,10)
fitness_vs_covs_lm_objects = list()
for (fitness_score_ind in fitness_scores_inds){
  y = fitness_scores[,fitness_score_ind]
  currname = colnames(fitness_scores)[fitness_score_ind]
  if(is.element(currname,set=names(fitness_vs_covs_lm_objects))){next}
  # All covariates, without the additionally excluded
  lm1 = get_lm_residuals(y,covariate_matrix,use_categorical=T,max_num_classes=20,
                         feature_is_numeric,additional_subcategories_to_exclude=potential_subcategories_to_exclude,
                         feature2subcategory=initial_cov_matrix_feature_subcategory)
  # # All covariates
  # lm2 = get_lm_residuals(y,covariate_matrix,use_categorical=T,max_num_classes=20,
  #                        feature_is_numeric,additional_subcategories_to_exclude=NULL)
  # 
  # plot(lm1[[1]]$residuals,lm2[[1]]$residuals)
  # par(mfrow=c(1,1),mar=c(4,4,4,4))
  plot(y[names(lm1[[1]]$residuals)],lm1[[1]]$residuals,ylab="Residual",xlab="Fitness score")
  lm1_summ = summary(lm1[[1]])
  sigs = lm1_summ$coefficients[,4]
  sort(sigs)[1:5]
  fitness_vs_covs_lm_objects[[currname]] = lm1
  break
}

gwas_cov_matrix = covariate_matrix
completion_status = pheno_data[,"Completion status of test.0.0"]
save(fitness_scores,gwas_cov_matrix,initial_cov_matrix,completion_status,fitness_vs_covs_lm_objects,file="UKBB_phenotypic_data_for_GWAS.RData")
#save(fitness_vs_covs_lm_objects,file="test_lm_size.RData")

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
is_higher_better = c(F,T)
names(is_higher_better) = names(fitness_vs_covs_lm_objects)
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
save(gwas_data,file="may31_2017_gwas_data_table.RData")
plot(gwas_data_residuals[,1],gwas_data_residuals[,2])

write.table(gwas_data,file = "may31_2017_gwas_data_table.txt",sep="\t",quote=F)

# Technical
# centre = "UK Biobank assessment centre.0.0"
# batch = get_regex_cols(all_features,"batch")
# get_regex_cols(all_features,"batch")
# batch_columns = get_regex_cols(all_features,"^genotype",ignore.case=T)
# apply(subject_cleaned_pheno_data[,batch_columns],2,table)
# length(get_regex_cols(all_features,"ECG"))

##########################################################

######## Analysis of the fitness scores ##########
library(corrplot)
par(mfrow=c(1,1))
corrplot(get_pairwise_corrs(fitness_scores[fully_comps,]))
corrplot(get_pairwise_corrs(fitness_scores[setdiff(rownames(fitness_scores),fully_comps),]))
par(mfrow=c(2,2),mar=c(2,4,2,2))
boxplot(get_list_by_values(subject_ols_preds_100,subject_technical_class,1,"No data"),main="HR_pred 100")
#boxplot(get_list_by_values(max_WD_pheno_data,subject_technical_class,1,"No data"), main="Max_WD")
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

# Chi square tests
Mc = initial_cov_matrix
Mf = fitness_scores
cut_size = 10
Mf = apply(Mf,2,cut_by_quantiles,nbreaks=cut_size)
for(j in 1:ncol(Mc)){
  if(feature_is_numeric[colnames(Mc)[j]]){
    v = as.numeric(as.character(Mc[,j]))
    v = cut_by_quantiles(v,nbreaks=min(cut_size,length(unique(v))))
    Mc[,j] = v
  }
}
chisq_pval_mat = matrix(1,nrow=ncol(Mf),ncol=ncol(Mc))
colnames(chisq_pval_mat) = colnames(Mc)
rownames(chisq_pval_mat) = colnames(Mf)
for(j1 in 1:ncol(Mf)){
  print(j1)
  for(j2 in 1:ncol(Mc)){
    v1 = Mf[,j1]
    v2 = Mc[,j2]
    inds = !is.na(v1) & !is.na(v2)
    tb = table(v1[inds],v2[inds])
    chisq_pval_mat[j1,j2] = chisq.test(tb)$p.value
  }
}
sort(chisq_pval_mat[2,])[1:50]
sort(chisq_pval_mat[10,])[1:20]

na_pval_mat = matrix(1,nrow=ncol(Mf),ncol=ncol(Mc))
colnames(na_pval_mat) = colnames(Mc)
rownames(na_pval_mat) = colnames(Mf)
for(j1 in 1:ncol(Mf)){
  print(j1)
  for(j2 in 1:ncol(Mc)){
    v1 = is.na(Mf[,j1])
    v2 = Mc[,j2]
    inds = !is.na(v2)
    tb = table(v1[inds],v2[inds])
    na_pval_mat[j1,j2] = chisq.test(tb)$p.value
  }
}
sort(na_pval_mat[1,])[1:20]

save(na_pval_mat,chisq_pval_mat,file = "fitness_scores_vs_features_chisq_analysis.RData")

# another analysis type - centered for a fitness score
fitness_score_inds = c(2,9,10)
cut_size = 5
Mc = initial_cov_matrix
for(j in 1:ncol(Mc)){
  if(feature_is_numeric[colnames(Mc)[j]]){
    v = as.numeric(as.character(Mc[,j]))
    v = cut_by_quantiles(v,nbreaks=min(cut_size,length(unique(v))))
    Mc[,j] = v
  }
}

library(entropy)
covariance_correlation_summary_tables = list()
for(fit_ind in fitness_score_inds){
  fv = fitness_scores[,fit_ind]
  fv_disc = cut_by_quantiles(fv,nbreaks = cut_size)
  fv_is_na = is.na(fv)
  summary_table = c()
  for(j in 1:ncol(Mc)){
    covv = initial_cov_matrix[,j]
    covv_disc = Mc[,j]
    covv_na = is.na(covv)
    cov_name = colnames(Mc)[j]
    cov_category = initial_cov_matrix_feature_subcategory[cov_name]
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
}
corrected_names = c("HR_exercise_pred_100","HR_recovery_ratio_60","max_WD")
names(covariance_correlation_summary_tables) = corrected_names
for(nn in corrected_names){
  colnames(covariance_correlation_summary_tables[[nn]]) = c("Feature","Category",
                                                            "MI-discretized","ChisqP-discretized",
                                                            "MI-NA fitness","ChisqP-NA fitness",
                                                            "Spearman rho","Spearman rho p",
                                                            "MI-NAs", "ChisqP-NAs")
  write.table(covariance_correlation_summary_tables[[nn]],file=paste(nn,"_cov_analysis.txt",sep=""),sep="\t",quote=F,row.names=F)
}
save(covariance_correlation_summary_tables,file="covariance_correlation_summary_tables.RData")

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
