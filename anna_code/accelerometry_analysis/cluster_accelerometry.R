rm(list=ls())
data=read.table("all.features.na.omit.nooutliers.txt",header=TRUE,sep='\t',stringsAsFactors=FALSE)
rownames(data)=as.character(data$Subject) 
data$Subject=NULL

#k-means with high value of k 
k=100 
fit=kmeans(data,k)
centroids=fit$centers
#hierarchical clustering on the centroids 
d <- dist(centroids, method = "euclidean") # distance matrix
hc <- hclust(d, method="average")
plot(hc) # display dendogram
#groups <- cutree(fit, k=5) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters
#rect.hclust(fit, k=5, border="red") 