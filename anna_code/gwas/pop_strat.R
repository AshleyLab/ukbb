rm(list=ls())
library(ggplot2)
data=read.table("pca_results_genome.filtered.eigenvec",header=TRUE,sep='\t')
#color-code by self-reported ancestry 
ancestry=read.table("ethnicity_aggregate_phenotypes.categorical.converted.txt",header=TRUE,sep='\t')
merged=merge(data,ancestry,by="FID")
#plot the prinicipal components! 
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
colors=getPalette(22)
ggplot(data=merged,aes(x=merged$PC1,
                       y=merged$PC2,
                       color=merged$Ethnicity))+
  geom_point()+
  xlab("PC1")+
  ylab("PC2")+
  ggtitle("Batch 1 Pop Strat")+
  scale_color_manual(values = colors)
  theme_bw()


  
ggplot(data=merged,aes(x=merged$PC2,
                         y=merged$PC3,
                         color=merged$Ethnicity))+
    geom_point()+
    xlab("PC2")+
    ylab("PC3")+
    ggtitle("Batch 1 Pop Strat")+
    scale_color_manual(values = colors)
  theme_bw()
  
ggplot(data=merged,aes(x=merged$PC1,
                         y=merged$PC3,
                         color=merged$Ethnicity))+
    geom_point()+
    xlab("PC1")+
    ylab("PC3")+
    ggtitle("Batch 1 Pop Strat")+
    scale_color_manual(values = colors)
    theme_bw()

tokeep_FID=merged$FID[merged$PC1>-0.002]
tokeep_IID=merged$IID.x[merged$PC1 > -0.002]
tokeep=data.frame(tokeep_FID,tokeep_IID)
write.table(tokeep,file="pop_strat_filtered.txt",quote=FALSE,row.names = FALSE,col.names = FALSE,sep='\t')  
