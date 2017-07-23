setwd("/Users/David/Desktop/ukbb/tests")
xx = get(load("may31_2017_gwas_data_table.RData"))
yy = read.delim("Predicted_HR_at_WD=100_conservative_euro.txt",row.names=1,header=T)
inter = intersect(rownames(yy),rownames(xx))
par(mfrow=c(1,2))
plot(x=xx[inter,1],y=yy[inter,2],xlab="May",ylab="July",main="May vs. July raw")
#plot(x=xx[inter,1],y=yy[inter,2],xlab="May",ylab="July",xlim=c(-50,100),ylim=c(-50,100),main="May vs. July raw")
xx = xx[!is.na(xx[,1]),]
yy = yy[!is.na(yy[,2]),]
in_may_not_july = setdiff(rownames(xx),rownames(yy))
in_july_not_may = setdiff(rownames(yy),rownames(xx))
in_both = intersect(rownames(yy),rownames(xx))
# some tests before quantile
may_scores_in_both = xx[in_both,1]
july_scores_in_both = yy[in_both,2]

comparison_table = cbind(may_scores_in_both,july_scores_in_both)
write.table(comparison_table,file="comparison_table.txt",sep="\t",quote=F,col.names = T,row.names = T)

cor(may_scores_in_both,july_scores_in_both,method="spearman")
cor(may_scores_in_both,july_scores_in_both)
may_scores_in_may_not_july = xx[in_may_not_july,1]
wilcox.test(may_scores_in_both,may_scores_in_may_not_july)
july_scores_in_july_not_may = yy[in_july_not_may,2]
wilcox.test(july_scores_in_both,july_scores_in_july_not_may)
# some tests after quantile
library(preprocessCore)
xx[,1] = simple_normalization(xx[,1])
yy[,2] = simple_normalization(yy[,2])
may_scores_in_both = xx[in_both,1]
july_scores_in_both = yy[in_both,2]
plot(x=may_scores_in_both,y=july_scores_in_both,xlab="May",ylab="July",main="May vs. July after normalization")
wilcox.test(may_scores_in_both,july_scores_in_both,paried=T)
cor(may_scores_in_both,july_scores_in_both,method="spearman")
cor(may_scores_in_both,july_scores_in_both)
# may_scores_in_may_not_july = xx[in_may_not_july,1]
# wilcox.test(may_scores_in_both,may_scores_in_may_not_july)
# july_scores_in_july_not_may = yy[in_july_not_may,2]
# wilcox.test(july_scores_in_both,july_scores_in_july_not_may)
# Set overlaps
library(VennDiagram)
area1 = nrow(xx)
area2 = nrow(yy)
cross = length(in_both)

simple_normalization<-function(x){
  r = rank(x)
  r = r/(max(r)+1)
  z = qnorm(r)
  return(z)
}




