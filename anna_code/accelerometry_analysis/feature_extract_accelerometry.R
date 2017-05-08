library(wavelets) 
source('extract_features.R')
source('helpers.R')

args <- commandArgs(trailingOnly = TRUE)
datafiles=args[1]
datafiles<-read.table(datafiles,header=FALSE,sep=',')$V1
for(datafile in datafiles){
	     parts=unlist(strsplit(datafile,'/'))
	     last_part=parts[length(parts)]
	     subject=unlist(strsplit(last_part,'_'))[1]
	     data=na.omit(as.numeric(as.matrix(read.table(datafile,header=FALSE,sep=','))))
	     tryCatch({
		result_dwt<-c(subject,unlist(dwt_transform_features(data,subject,13)))
		result_dwt=as.matrix(t(result_dwt))
		write.table(result_dwt,file=paste(subject,"features.csv",sep='.'),quote=FALSE,sep='\t',row.names=FALSE,col.names=FALSE)
		},
		error=function(e){})
		}