data=read.table("afib.cohort.txt",header=FALSE) 
n_train=0.8*nrow(data) 
train_indices=sample(nrow(data),n_train,replace=FALSE)
train_split=data[train_indices,]
test_split=data[-train_indices,]
write.table(train_split,"afib.cohort.training.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(test_split,"afib.cohort.test.txt",quote=FALSE,row.names=FALSE,col.names=FALSE) 

