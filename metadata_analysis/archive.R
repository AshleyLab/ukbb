# This is an archive with old commented out code
# ##### Get the r2 values after taking a window of linear WD increases #######
# get_window_by_constraints<-function(x,y,max_dev=2){
#   n = length(x)
#   window_mat = matrix(F,nrow=n,ncol=n)
#   best_window = c(1,1)
#   for(i in 1:(n-1)){
#     for(j in (i+1):n){
#       a = x[i:j];b=y[i:j]
#       res = abs(lm(b~a)$residuals)
#       if(all(res<max_dev)){
#         window_mat[i,j]=T
#         if(j-i > (best_window[2]-best_window[1])){
#           best_window = c(i,j)
#         }
#       }
#     }
#   }
#   return(best_window)
# }
# 
# r2_values_in_subwindow = c()
# for(i in names(lm_inputs)){
#   if(is.null(lm_inputs[[i]])){next}
#   WD = lm_inputs[[i]]$WD; HR = lm_inputs[[i]]$HR; TP = lm_inputs[[i]]$TP
#   if (length (WD) < 3 ){next}
#   orig_size = length(WD)
#   w_inds = get_window_by_constraints(TP,WD)
#   WD = WD[w_inds[1]:w_inds[2]];TP = TP[w_inds[1]:w_inds[2]];HR = HR[w_inds[1]:w_inds[2]]
#   new_size = length(WD)
#   if (length (WD) < 2 ){next}
#   obj1 = lm(HR~WD); obj2 = lm(HR~TP); obj3 = lm(HR~TP+WD)
#   r2_values_in_subwindow = rbind(r2_values_in_subwindow,
#                                  c(get_lm_r2(obj1),get_lm_r2(obj2),get_lm_r2(obj2),orig_size,new_size)
#   )
#   rownames(r2_values_in_subwindow)[nrow(r2_values_in_subwindow)] = i
#   if(as.numeric(i)%%100==0){print(i);print(w_inds)}
# }
# plot(r2_values_in_subwindow[,1],r2[rownames(r2_values_in_subwindow)])

# # A technical measure: correlation of time vs. WD
# wd_tp_pearson_cor = c();wd_tp_spearman_cor = c()
# for(i in names(lm_inputs)){
#   if(is.null(lm_inputs[[i]])){next}
#   WD = lm_inputs[[i]]$WD; HR = lm_inputs[[i]]$HR;TP = lm_inputs[[i]]$TP
#   if(length(WD) < 3){next}
#   currn = length(wd_tp_spearman_cor)
#   wd_tp_pearson_cor[currn+1] = cor(WD,TP)
#   names(wd_tp_pearson_cor)[currn+1] = i
#   wd_tp_spearman_cor[currn+1] = cor(WD,TP,method="spearman")
#   names(wd_tp_spearman_cor)[currn+1]=i
#   if(as.numeric(i)%%1000==0){print(i)}
# }
# 
# # A technical measure that quantifies dynamic deviation from the experiment design
# subj2wd_tp_experiment_score = c()
# num_points = 3;last_points = 3
# for(i in names(lm_inputs)){
#   if(is.null(lm_inputs[[i]])){next}
#   #WD = regr_inputs[[i]]$WD; HR = regr_inputs[[i]]$HR;TP = regr_inputs[[i]]$TP
#   WD = lm_inputs[[i]]$WD; HR = lm_inputs[[i]]$HR;TP = lm_inputs[[i]]$TP
#   if(length(WD)<6){subj2wd_tp_experiment_score[i]=0;next}
#   tr_WD = WD[1:3];tr_TP = TP[1:3]
#   te_WD = WD[-c(1:3)]; te_TP = TP[-c(1:3)]
#   obj = lm(tr_WD~tr_TP)
#   preds = predict(obj,newdata=data.frame(tr_TP = te_TP))
#   devs = abs(preds - te_WD)
#   subj2wd_tp_experiment_score[i] = mean(devs[(length(devs)-last_points+1):(length(devs))])/mean(te_WD)
# }


# # The set of subject scores
# #wd_tp_cor = wd_tp_pearson_cor
# #linear_dev = subj2wd_tp_experiment_score
# m = cbind(r2[fully_comps],autoc[fully_comps],mean_HRs[fully_comps],sd_HRs[fully_comps])
# cor(m,method="spearman")
# 
# par(mfrow=c(1,1),mar=c(5,5,5,5)) 
# tble = c(table(r2[fully_comps]>0.5,autoc[fully_comps]>0.5))
# names(tble) = c("R^2<=0.5,Autoc<=0.5","R^2>0.5,Autoc<=0.5","R^2<=0.5,Autoc>0.5","R^2>0.5,Autoc>0.5")
# names(tble) = paste(names(tble),tble,sep="-")
# pie(tble)

# # Is the r2 correlated with any of the scores
# y1 = m[,1]
# x1 = m[,-c(1,2)]
# x1 = cbind(x1,x1[,1]/x1[,2])
# summary(lm(y1~log(x1)))




# hist(wd_tp_pearson_cor^2,breaks=200,xlim=c(0.97,1),main="WD vs. Time R2")
# subj2wd_vs_tp_class = rep("A",length(wd_tp_pearson_cor));names(subj2wd_vs_tp_class) = names(wd_tp_pearson_cor)
# subj2wd_vs_tp_class[wd_tp_pearson_cor^2 < 0.98] = "B"
# boxplot(get_list_by_values(r2,subj2wd_vs_tp_class[names(r2)]),ylab="R^2",las=2)

# # Look at some specific cases
# subjects = intersect(rownames(autocorrs),names(sizes))
# subject2class = subject_technical_class[subjects]
# fully_comps = intersect(subjects,names(which(subject2class=="C:Large")))
# 
# # Example 0: high r2, high autocorr, fully completed
# # Possible explanation: outliers or change in WD increase rate (e.g., two slopes)
# inds = which(autocorrs[fully_comps,1]>0.5 & r2[fully_comps]>0.5)
# print(length(inds))
# i = names(inds)[10000]
# HR = lm_inputs[[i]]$HR;TP = lm_inputs[[i]]$TP;WD = lm_inputs[[i]]$WD
# par(mfrow=c(1,2),mar = c(4,4,4,2))
# plot(TP,HR);plot(TP,WD)
# 
# # Example 1.2: low r2, high autocorr, fully completed
# # Possible explanation: outliers or change in WD increase rate (e.g., two slopes)
# inds = which(autocorrs[fully_comps,1]>0.5 & r2[fully_comps]<0.5)
# print(length(inds))
# i = names(inds)[13]
# HR = lm_inputs[[i]]$HR;TP = lm_inputs[[i]]$TP;WD = lm_inputs[[i]]$WD
# par(mfrow=c(1,2),mar = c(4,4,4,2))
# plot(TP,HR);plot(TP,WD)
# w_inds = get_window_by_constraints(TP,WD)
# HR = HR[w_inds[1]:w_inds[2]];WD = WD[w_inds[1]:w_inds[2]];TP=TP[w_inds[1]:w_inds[2]]
# summary(lm(HR~WD))
# acf(HR)
# 
# # Example 2: low r2, low autocorr, fully completed
# inds = which(autocorrs[fully_comps,1]<0.5 & r2[fully_comps]<0.2)
# print(length(inds))
# i = names(inds)[2]
# HR = lm_inputs[[i]]$HR;TP = lm_inputs[[i]]$TP;WD = lm_inputs[[i]]$WD
# par(mfrow=c(1,2),mar = c(4,4,4,2))
# plot(TP,HR);plot(TP,WD)
# w_inds = get_window_by_constraints(TP,WD)
# HR = HR[w_inds[1]:w_inds[2]];WD = WD[w_inds[1]:w_inds[2]]
# summary(lm(HR~WD))
# 
# # Example 3: high r2, low autocorr, fully completed
# inds = which(autocorrs[fully_comps,1]<0.1 & r2[fully_comps]>0.6)
# print(length(inds))
# i = names(inds)[1]
# HR = lm_inputs[[i]]$HR;TP = lm_inputs[[i]]$TP;WD = lm_inputs[[i]]$WD
# par(mfrow=c(1,2),mar = c(4,4,4,2))
# plot(TP,HR);plot(TP,WD)

# IN PROGRESS:
# Back-test analysis and RMSE analysis
# This should be used as an extra validation for the analysis
# We compare different methods in their ability to predict the last k points
# of an individual.
# Since we saw that the workloads increase at a constant rate then using the "next k" can be 
# used to learn the heart rate at a specific point and vice versa
# The compared methods
# 1. Simple linear regression of workloads vs. heart rate. 
#       Different HRs of the same workloads are not averaged
# 2. Simple linear regression of time point vs. heart rate (e.g., simple trend analysis). 
#       Different HRs of the same workloads are not averaged
# 3. Non-stationary time series analysis using ARMA
#       Different HRs of the same workloads are averaged
# This is a wrapper that runs a regression analysis for a given subject
# The training data are all time points outside a user given window
# The time points for the tests end at window_start (set to the last one if -1)
# and contain the preceding k time points.
# The regression algorithm is run on the training and the results on the test window are
# the output.
# The function assumes that x and y has at least three instances
# run_last_points_regression_validation<-function(y,x,regression_algo,window_size=1,
#                                                 window_end=-1,remove_duplications=F,...){
#   if(window_end<0){window_end = length(y)}
#   if(is.null(dim(x))){x = matrix(x,ncol=1)}
#   if(is.null(colnames(x))){colnames(x) = paste("V",1:ncol(x),sep="")}
#   in_test = rep(F,length(y))
#   in_test[(window_end-window_size+1):window_end] = T
#   # Deal with duplications
#   if(remove_duplications){
#     in_test_y_values = unique(y[in_test])
#     in_test = is.element(y,set=in_test_y_values)
#   }
#   test_set_x = x[in_test,];test_set_y=y[in_test]
#   tr_set_x = x[!in_test,];tr_set_y=y[!in_test]
#   if(is.null(dim(test_set_x))){test_set_x = matrix(test_set_x,ncol=1)}
#   if(is.null(dim(tr_set_x))){tr_set_x = matrix(tr_set_x,ncol=1)}
#   colnames(test_set_x) = colnames(x)
#   colnames(tr_set_x) = colnames(x)
#   if(!is.null(names(y))){
#     rownames(x) = names(y)
#     rownames(test_set_x) = rownames(x)[in_test]
#     rownames(tr_set_x) = rownames(x)[!in_test]
#   }
#   regr_obj = regression_algo(tr_set_y,tr_set_x,...)
#   predictions = predict(regr_obj,data.frame(test_set_x))
#   real_values = test_set_y
#   return(cbind(predictions,real_values))
# }
# 
# # This is a wrapper that runs a time series modelling analysis for a given subject
# # The training data are all time points outside a user given window
# # The time points for the tests end at window_start (set to the last one if -1)
# # and contain the preceding k time points.
# # The time series algorithm is run on the training and the forecasting 
# # on the test window are the output.
# # The function assumes that y has at least three instances
# # And there is no need to exlude duplications in y (i.e., we average HRs with the same WD)
# # If x is NULL: run ARIMA on y, otherwise use VARMA
# run_last_points_arima_validation<-function(y,x=NULL,arima_func=arima0,window_size=1,window_end=-1,...){
#   if(window_end<0){window_end = length(y)}
#   in_test = rep(F,length(y))
#   in_test[(window_end-window_size+1):window_end] = T
#   test_set_y=y[in_test];tr_set_y=y[!in_test]
#   if(is.null(x)){
#     tr_obj = arima_func(tr_set_y,...)
#   }
#   else{
#     tr_obj = arima_func(tr_set_y,x,...)
#   }
#   predictions = as.numeric(predict(tr_obj,length(test_set_y))$pred)
#   real_values = test_set_y
#   return(cbind(predictions,real_values))
# }
# get_preds_mse <- function(res){return(mean((res[,1]-res[,2])^2))}
# 
# # library(MTS)
# # run_VARX<-function(y,x,y_order=1,x_order=1,...){
# #   obj = list()
# #   obj[[1]] = VARX(y,y_order,x,x_order,...)
# #   class(obj) = "run_VARX"
# #   return(obj)
# # }
# # predict.run_VARX<-function(obj,x,with_se=F){
# #   preds = VARXpred2(obj[[1]],x)
# #   if(with_se){return(preds)}
# #   if(is.null(dim(preds))){return(preds[1])}
# #   return(preds[,1])
# # }
# # # Tests
# # WD = lm_inputs[[i]]$WD
# # HR = lm_inputs[[i]]$HR
# # TP = lm_inputs[[i]]$T
# # HR = as.numeric(tapply(HR,WD,mean))
# # TP = as.numeric(tapply(TP,WD,mean))
# # WD = as.numeric(tapply(WD,WD,mean))
# # lm_pred = run_last_points_regression_validation(HR,WD,run_lm,5)
# # ar_pred = run_last_points_arima_validation(HR,window_size=5,order = c(5,0,5))
# # get_preds_mse(lm_pred)
# # get_preds_mse(ar_pred)
# # # ARIMA models
# # library(stats)
# # ar_model = arima0(HR,order = c(5,0,0))
# # ma_model = arima0(HR,order = c(0,0,5))
# # predict(ar_model,5)$pred
# # varx_obj = run_VARX(HR,WD)
# # VARXpred(varx_obj[[1]],newxt = 5)
# # predict(varx_obj,114)
# # varx_pred = run_last_points_regression_validation(HR,WD,run_VARX,2)
# 
# # Back test analysis (cross validation like analysis)
# window_size = 2
# algo2pred_vs_real = list()
# subject_MSEs = c()
# for(i in 1:length(lm_inputs)){
#   if(is.null(lm_inputs[[i]])){next}
#   print(i)
#   WD = lm_inputs[[i]]$WD
#   HR = lm_inputs[[i]]$HR
#   TP = lm_inputs[[i]]$TP
#   if(length(WD)<(2*window_size)){next}
#   res_lst = list()
#   res_lst[["lm_1"]] = run_last_points_regression_validation(HR,WD,run_lm,window_size)
#   res_lst[["lm_2"]] = run_last_points_regression_validation(HR,TP,run_lm,window_size)
#   #try({res_lst[["arima_5_0_0"]] = run_last_points_arima_validation(HR,window_size=window_size,order=c(2,0,0))})
#   for(nn in names(res_lst)){
#     rownames(res_lst[[nn]]) = paste(names(lm_inputs)[i],1:nrow(res_lst[[nn]]),sep=";")
#     algo2pred_vs_real[[nn]] = rbind(algo2pred_vs_real[[nn]],res_lst[[nn]])
#   }
#   subject_MSEs=rbind(subject_MSEs,sapply(res_lst,get_preds_mse))
#   rownames(subject_MSEs)[nrow(subject_MSEs)] = names(lm_inputs)[i]
# }
# colMeans(sqrt(subject_MSEs),na.rm=T)
# for(j in 1:length(algo2pred_vs_real)){
#   m = algo2pred_vs_real[[j]]
#   m = m[grepl(";1",rownames(m)),]
#   print(cor(m))
# }
# subject_rmse1 = sqrt(subject_MSEs)[,1]
# #############################################################
# ## Look at some examples of subjects with very high RMSE
# # The subject with the worst RMSE
# table(subject_rmse1>100)
# subj_example1 = names(which(subject_rmse1>100)[1])
# subject_r_squared[subj_example1]
# subj_example1_HR = lm_inputs[[subj_example1]]$HR
# subj_example1_WD = lm_inputs[[subj_example1]]$WD
# subj_example1_TP = lm_inputs[[subj_example1]]$T
# plot(subj_example1_TP,subj_example1_HR)
# plot(subj_example1_WD)
# summary(lm(subj_example1_HR~subj_example1_WD))
# subject_ols_preds[subj_example1]
# subject_MSEs[subj_example1,]^0.5
# # Bad RMSE, high R2
# subj_example2 = names(which(subject_r_squared<0.2)[2])
# subject_r_squared[subj_example2]
# subj_example2_HR = lm_inputs[[subj_example2]]$HR
# subj_example2_WD = lm_inputs[[subj_example2]]$WD
# subj_example2_TP = lm_inputs[[subj_example2]]$T
# plot(lm(subj_example2_HR~subj_example2_TP))
# plot(subj_example2_WD)
# plot(subj_example2_WD,subj_example2_HR)
# plot(y=subj_example2_HR,subj_example2_TP)
# summary(lm(subj_example2_HR~subj_example2_WD))
# subject_ols_preds[subj_example2]
# subject_MSEs[subj_example1,]^0.5
# 
# ##### TODO: check these subjects: they have workloads but they are not monotone
# table(not_mono_with_wlds)
# rand_ind = sample(which(not_mono_with_wlds))[1]
# plot(new_time_order[rand_ind,],regex_matrices[["workload"]][rand_ind,])
# plot(regex_matrices[["heartrate"]][rand_ind,],regex_matrices[["workload"]][rand_ind,1:112])
# tmp_x = as.numeric(regex_matrices[["heartrate"]][rand_ind,])
# tmp_y = as.numeric(regex_matrices[["workload"]][rand_ind,1:112])
# inds = !is.na(tmp_x) & !is.na(tmp_y)
# tmp_x = tmp_x[inds];tmp_y=tmp_y[inds]
# tmp_x = tmp_x + rnorm(length(tmp_x))/10
# tmp_y = tmp_y + rnorm(length(tmp_y))/10
# cor.test(tmp_x,tmp_y,method="kendall")
# cor(tmp_x,tmp_y,method="pearson")
# summary(lm(tmp_y~tmp_x))
# plot(tmp_x,tmp_y)
# lm_obj = lm(tmp_y~tmp_x)
# res = lm_obj$residuals
# plot(res)
# 
# # Test time series multivariate analysis
# library(MTS)
# var_model = VAR(cbind(tmp_x,tmp_y),p=2)
# plot(var_model$residuals[,2])
# plot(var_model$residuals)
# median(abs(var_model$residuals[,2]))
# median(abs(lm_obj$residuals))
# test_last_points<-function()
# ########################################


# Additional regression analyses, some of which are time-series models
# library(MTS)
# VAR_models = list()
# for (j in 1:length(tps_vs_wkls)){
#   print(j)
#   if(is.null(lm_inputs[[j]])){next}
#   x = lm_inputs[[j]]$x
#   y = lm_inputs[[j]]$y
#   inds = !is.na(x) & !is.na(y)
#   x = x[inds];y=y[inds]
#   try({
#     var_model = VAR(cbind(y,x),p=2)
#     VAR_models[[j]] = var_model
#   })
# }
# get_r_squared <-function(obj){
#   if (is.null(obj)){return (-1)}
#   return(summary(obj)[[8]])
# }
# get_residual_stats<-function(obj,func=mean){
#   if (is.null(obj)){return (-1)}
#   return (func(obj$residuals))
# }
# 
# get_VAR_residual_stats<-function(obj,func=median,col=1){
#   if (is.null(obj)){return (-1)}
#   return (func(obj$residuals[,col]))
#}
# var_med_abs_res = sapply(VAR_models,get_VAR_residual_stats,func = function(x)median(abs(x)))
# lm_med_abs_res = sapply(lm_objs,get_residual_stats,func = function(x)median(abs(x)))
# plot(var_med_abs_res,lm_med_abs_res);abline(0,1)
# diffs = abs(lm_med_abs_res - var_med_abs_res)
# max_ind = sample(which(diffs>15))[1]
# max_ind_x = lm_inputs[[max_ind]]$x
# max_ind_y = lm_inputs[[max_ind]]$y
# par(mfrow=c(2,2))
# plot(max_ind_y)
# plot(max_ind_x)
# plot(max_ind_x,max_ind_y)
# VARpred(VAR_models[[max_ind]],1)
# plot(lm_objs[[max_ind]]$residuals)
# 
# max_ind_TPs = 1:length(max_ind_x)
# new_lm = lm (max_ind_y~max_ind_TPs+max_ind_x)
# plot(new_lm$residuals,lm_objs[[max_ind]]$residuals)

# # OLD code - we looked at some classes manually in order to better understand the database
# #  A series of covariates by name type
# pheno_data_cov_type_names = list()
# pheno_data_cov_submatrices = list()
# pheno_data_cov_type_names[["genetic_pca_mat"]] = get_regex_cols(all_features,"genetic principal",ignore.case=T)
# 
# # Test related
# pheno_data_cov_type_names[["exercise_test_features"]] = c("Duration of fitness test.0.0","Completion status of test.0.0",
#                            "Test completion status.0.0",
#                            "Maximum workload during fitness test.0.0",
#                            "Maximum heart rate during fitness test.0.0",
#                            "Duration of fitness test.0.0",
#                            "Target heart rate achieved.0.0",
#                            "Chest pain felt during physical activity.0.0")
# #pheno_data_cov_type_names[["diet_features"]] = c("Duration of fitness test.0.0","Completion status of test.0.0")
# # Result ranking?
# pheno_data_cov_type_names[["result_ranking_features"]] = get_regex_cols(all_features,"ranking",ignore.case=T)
# # Blood cell tests
# blood_cell_features = get_regex_cols(all_features,"blood cell",ignore.case=T)
# blood_cell_features = c(blood_cell_features,get_regex_cols(all_features,"lympho",ignore.case=T))
# blood_cell_features = c(blood_cell_features,get_regex_cols(all_features,"monocy",ignore.case=T))
# blood_cell_features = c(blood_cell_features,get_regex_cols(all_features,"phil",ignore.case=T))
# blood_cell_features = c(blood_cell_features,get_regex_cols(all_features,"haemo",ignore.case=T))
# blood_cell_features = c(blood_cell_features,get_regex_cols(all_features,"reticulocyte",ignore.case=T))
# blood_cell_features = c(blood_cell_features,get_regex_cols(all_features,"platelet",ignore.case=T))
# blood_cell_features = c(blood_cell_features,get_regex_cols(all_features,"sphered cell",ignore.case=T))
# blood_cell_features = c(blood_cell_features,get_regex_cols(all_features,"Haematocrit",ignore.case=T))
# blood_cell_features = c(blood_cell_features,get_regex_cols(all_features,"corpuscular",ignore.case=T))
# blood_cell_features = c(blood_cell_features,get_regex_cols(all_features,"sphered cell",ignore.case=T))
# blood_cell_features = unique(blood_cell_features)
# pheno_data_cov_type_names[["blood_cell_features"]] = blood_cell_features
# # Fat measure features
# fat_features = get_regex_cols(all_features," fat ",ignore.case=T)
# fat_features = union(fat_features,get_regex_cols(all_features," fat-",ignore.case=T))
# pheno_data_cov_type_names[["fat_features"]] = fat_features
# pheno_data_cov_type_names[["device_features"]] = get_regex_cols(all_features,"device",ignore.case=T)
# pheno_data_cov_type_names[["predicted_scores"]] = get_regex_cols(all_features,"predicted",ignore.case=T)
# pheno_data_cov_type_names[["impedance"]] = get_regex_cols(all_features,"impedance",ignore.case=T)
# pheno_data_cov_type_names[["mass"]] = get_regex_cols(all_features,"mass",ignore.case=T)
# cognitive_features = get_regex_cols(all_features,"snap-button",ignore.case=T)
# cognitive_features = union(cognitive_features,get_regex_cols(all_features,"matches in round",ignore.case=T))
# cognitive_features = union(cognitive_features,get_regex_cols(all_features,"Time to complete round",ignore.case=T))
# pheno_data_cov_type_names[["cognitive_features"]] = cognitive_features
# # Get the matrices and remove them from the cov matrix
# for(mat_name in names(pheno_data_cov_type_names)){
#   m = pheno_data[,pheno_data_cov_type_names[[mat_name]]]
#   covariate_matrix = covariate_matrix[,setdiff(colnames(covariate_matrix),colnames(m))]
#   pheno_data_cov_submatrices[[mat_name]] = m
# }
# dim(covariate_matrix)
# colnames(covariate_matrix)
# Covariates to exclude and features that we must have

# covs_to_exclude = c(
#   "Spirometry method.0.0",
#   "Date of attending assessment centre.0.0",
#   get_regex_cols(all_features,"predicted",ignore.case=T)
# )
# must_have_covs = c(
#   "Sex.0.0",
#   "UK Biobank assessment centre.0.0",
#   "Age when attended assessment centre.0.0",
#   "Genotype measurement batch.0.0",
#   "Body mass index (BMI).0.0",
#   "Number of self-reported cancers.0.0",
#   "Number of self-reported non-cancer illnesses.0.0",
#   "Number of operations, self-reported.0.0",
#   "Current tobacco smoking.0.0",
#   "Weight.0.0",
#   "Basal metabolic rate.0.0",
#   "Ethnic background.0.0",
#   "Ever smoked.0.0",
#   "Alcohol drinker status.0.0",
#   "Smoking status.0.0",
#   "Types of physical activity in last 4 weeks.0.0",
#   "Blood clot, DVT, bronchitis, emphysema, asthma, rhinitis, eczema, allergy diagnosed by doctor.0.0",
#   "Vascular/heart problems diagnosed by doctor.0.0",
#   "Mouth/teeth dental problems.0.0",
#   "Plays computer games.0.0",
#   "Long-standing illness, disability or infirmity.0.0",
#   "Overall health rating.0.0",
#   "Skin colour.0.0",
#   "Exposure to tobacco smoke at home.0.0",
#   "Sleep duration.0.0",
#   "Number of days/week of vigorous physical activity 10+ minutes.0.0",
#   "Duration of walks.0.0",
#   "Average total household income before tax.0.0",
#   "Townsend deprivation index at recruitment.0.0"
# )
# # Comment - "Current tobacco smoking.0.0" is not ordered properly, need to correct

# # OLD: other covariate analyses
# # Correlation between fitness scores and the covariates
# all(rownames(covariate_matrix) == rownames(fitness_scores))
# cov_fitness_corrs = c()
# for(i in 1:ncol(initial_cov_matrix)){
#   x = as.numeric(as.character(initial_cov_matrix[,i]))
#   if(all(is.na(x))){next}
#   for(j in 1:ncol(fitness_scores)){
#     y = fitness_scores[,j]
#     inds = !is.na(x) & !is.na(y) & !is.nan(x) & !is.nan(y)
#     currname = paste(colnames(initial_cov_matrix)[i],colnames(fitness_scores)[j],sep=";")
#     cov_fitness_corrs[currname] = cor(x[inds],y[inds])
#     print(cov_fitness_corrs[currname] )
#   }
# }
# 
# # NA analysis
# discrete_fitness_scores = apply(fitness_scores[,6:7],2,cut,breaks=10)
# 
# all(rownames(covariate_matrix) == rownames(fitness_scores))
# cov_fitness_corrs = c()
# for(i in 1:ncol(initial_cov_matrix)){
#   x = as.numeric(as.character(initial_cov_matrix[,i]))
#   if(all(is.na(x))){next}
#   for(j in 1:ncol(fitness_scores)){
#     y = fitness_scores[,j]
#     inds = !is.na(x) & !is.na(y) & !is.nan(x) & !is.nan(y)
#     currname = paste(colnames(initial_cov_matrix)[i],colnames(fitness_scores)[j],sep=";")
#     cov_fitness_corrs[currname] = cor(x[inds],y[inds])
#     print(cov_fitness_corrs[currname] )
#   }
# }
# 
# table(is.na(cvs))
# sort(cov_fitness_corrs,decreasing = T)[1:50]
# table(pheno_data[,get_regex_cols(all_features,"cancer",ignore.case=T)[1]])
# cov_lm_objs = list()
# for(nn in colnames(covariate_matrix)){
#   try({
#     x_nn = covariate_matrix[,nn]
#     print (paste(nn,class(x_nn),mode(x_nn),length(unique(x_nn))))
#     cov_lm_objs[[nn]] = 
#       lm(subject_maximal_workload_y~.,data=data.frame(cbind(subject_maximal_workload_y,x_nn)))
#   })
# }


#### Old comments
# TODO: use the data above to calculate the slopes
# For each subject take the time points from the exercise phase.
# In addition, exclude the first time points for which the workload 
# was at "base" levels.
# Then, for the selected time points obtain thw workloads and the heart
# rates and rerun the old flow (to get a VO2max score).
##################

