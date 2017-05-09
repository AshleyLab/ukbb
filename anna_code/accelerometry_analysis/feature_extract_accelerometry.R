rm(list=ls())
library(wavelets) 
library(data.table)
source('extract_features.R')
source('helpers.R') 

#For parallel execution on a cluster! These are start and end indices for subject files to parse 

#args <- commandArgs(trailingOnly = TRUE)
#datafile<-args[1] 
#subject<-args[2] 
datafile="1129255_90004_0_0.csv"
subject="1129255"
#############################################################################################
#CREATE A DATA FRAME FOR FEATURE STORAGE
data=as.matrix(read.table(datafile,header=FALSE,sep=','))
result_dwt<-dwt_transform_features(data,subject,13)
write.csv(t(a),file="test.txt",quote=FALSE,sep="\t",append=TRUE,row.names=FALSE,col.names=FALSE)