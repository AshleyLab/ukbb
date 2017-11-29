# Functional analysis of the snp results file
try({setwd("/Users/david/Desktop/ukbb/")})

library(Vennerable);library(ggplot2)

# ############## MAF by 0.001 #############
# # Paths to the FUMA results for interpretation
# fuma_pooled_analysis_main_path = list(
#   exercise= "gwas_interpretation/fuma/filtered_maf_0.01_1000genomes/pooled_p1_5e8_p2_1e4_ld_0.6_without_maf/",
#   recovery = "gwas_interpretation/fuma/filtered_maf_0.01_1000genomes/recovery_pooled_p1_5e8_p2_1e4_ld_0.6_without_maf/"
# )
# fuma_pooled_analysis_maffilter_path = list(
#   exercise = "gwas_interpretation/fuma/filtered_maf_0.01_1000genomes/pooled_p1_5e8_p2_1e4_ld_0.6_with_maf/",
#   recovery = "gwas_interpretation/fuma/filtered_maf_0.01_1000genomes/recovery_pooled_p1_5e8_p2_1e4_ld_0.6_with_maf/"
# )
# saturation_analysis_paths = list(
#   exercise = list(
#     "33" = "gwas_interpretation/fuma/filtered_maf_0.01_1000genomes/33per_p1_5e8_p2_1e4_ld_0.6_without_maf/",
#     "67" = "gwas_interpretation/fuma/filtered_maf_0.01_1000genomes/67per_p1_5e8_p2_1e4_ld_0.6_without_maf/",
#     "75" = "gwas_interpretation/fuma/filtered_maf_0.01_1000genomes/75per_p1_5e8_p2_1e4_ld_0.6_without_maf/",
#     "80" = "gwas_interpretation/fuma/filtered_maf_0.01_1000genomes/80per_p1_5e8_p2_1e4_ld_0.6_without_maf/",
#     "85" = "gwas_interpretation/fuma/filtered_maf_0.01_1000genomes/85per_p1_5e8_p2_1e4_ld_0.6_without_maf/",
#     "90" = "gwas_interpretation/fuma/filtered_maf_0.01_1000genomes/90per_p1_5e8_p2_1e4_ld_0.6_without_maf/",
#     "100" = "gwas_interpretation/fuma/filtered_maf_0.01_1000genomes/pooled_p1_5e8_p2_1e4_ld_0.6_without_maf/"
#   )
# )
# replication_analysis_files = list(
#   exercise = list(
#     "rep" = "gwas_interpretation/fuma_in_files/filtered_maf_0.001/filtered_pval_0.01/hr_rep_data.txt",
#     "67" = "gwas_interpretation/fuma/filtered_maf_0.01_1000genomes/67per_p1_5e8_p2_1e4_ld_0.6_without_maf/"
#   ),
#   recovery = list(
#     "rep" = "gwas_interpretation/fuma_in_files/filtered_maf_0.001/filtered_pval_0.01/recovery_rep_data.txt",
#     "67" = "gwas_interpretation/fuma/filtered_maf_0.01_1000genomes/recovery_67per_p1_5e8_p2_1e4_ld_0.6_without_maf/"
#   )
# )

############## MAF by 0.01 #############
# Paths to the FUMA results for interpretation
fuma_pooled_analysis_main_path = list(
  exercise= "gwas_interpretation/fuma/filtered_maf_0.01/pooled_exercise_hr/",
  recovery = "gwas_interpretation/fuma/filtered_maf_0.01/pooled_recovery/",
  rhr = "gwas_interpretation/fuma/filtered_maf_0.01/pooled_rhr_nonfitness/"
)
fuma_pooled_analysis_maffilter_path = list(
  exercise= "gwas_interpretation/fuma/filtered_maf_0.01/pooled_exercise_hr/",
  recovery = "gwas_interpretation/fuma/filtered_maf_0.01/pooled_recovery/",
  rhr = "gwas_interpretation/fuma/filtered_maf_0.01/pooled_rhr_nonfitness/"
)
saturation_analysis_paths = list(
  exercise = list(
    "33" = "gwas_interpretation/fuma/filtered_maf_0.01/33per_exercise_hr/",
    "67" = "gwas_interpretation/fuma/filtered_maf_0.01/67per_exercise_hr/",
    "75" = "gwas_interpretation/fuma/filtered_maf_0.01/75per_exercise_hr/",
    "80" = "gwas_interpretation/fuma/filtered_maf_0.01/80per_exercise_hr/",
    "85" = "gwas_interpretation/fuma/filtered_maf_0.01/85per_exercise_hr/",
    "90" = "gwas_interpretation/fuma/filtered_maf_0.01/90per_exercise_hr/",
    "100" = "gwas_interpretation/fuma/filtered_maf_0.01/pooled_exercise_hr/"
  )
)
replication_analysis_files = list(
  exercise = list(
    "rep" = "gwas_interpretation/fuma_in_files/filtered_maf_0.01/filtered_pval_0.01/hr_rep_data.txt",
    "67" = "gwas_interpretation/fuma/filtered_maf_0.01/67per_exercise_hr/"
  ),
  recovery = list(
    "rep" = "gwas_interpretation/fuma_in_files/filtered_maf_0.01/filtered_pval_0.01/recovery_rep_data.txt",
    "67" = "gwas_interpretation/fuma/filtered_maf_0.01/67per_recovery/"
  )
)

# MAF data
mafs = read.delim("gwas_interpretation/euro_freq.txt",row.names = NULL,header=F,stringsAsFactors=F,sep=" ")
snp2maf = as.numeric(mafs[,3])
names(snp2maf) = mafs[,2]
rm(mafs);gc()
snp2maf = pmin(snp2maf,1-snp2maf);gc()

##############################################
##############################################
##############################################
# Some functions for loading data from the FUMA output files and for replication analysis
get_fuma_gwascatalog_results<-function(path){
  res = read.delim(paste(path,"gwascatalog.txt",sep=''),stringsAsFactors = F)
  res = res[,c("GenomicLocus","IndSigSNP","snp","Trait","Region","ReportedGene","MappedGene","SNPs","P")]
  return(res)
}
get_fuma_lead_snp_data<-function(path){
  res = read.delim(paste(path,"leadSNPs.txt",sep=''),stringsAsFactors = F)
  rownames(res) = res$rsID
  return(res)
}
get_fuma_all_snp_data<-function(path){
  res = read.delim(paste(path,"snps.txt",sep=''),stringsAsFactors = F)
  rownames(res) = res$rsID
  return(res)
}
get_fuma_mapped_gene_data<-function(path){
  res = read.delim(paste(path,"genes.txt",sep=''),stringsAsFactors = F,row.names = 1)
  return(res)
}
get_fuma_gene_sets_results<-function(path){
  res = read.delim(paste(path,"GS.txt",sep=''),stringsAsFactors = F)
  return(res)
}
library(org.Hs.eg.db)
entrez2symbol = as.list(org.Hs.egSYMBOL)
symbol2entrez = as.list(org.Hs.egSYMBOL2EG)
ensembl2entrez = as.list(org.Hs.egENSEMBL2EG)
entrez2ensmbl = as.list(org.Hs.egENSEMBL)
ensembl2symbol = sapply(ensembl2entrez,function(x,y)unique(y[x]),y=entrez2symbol)
symbol2ensembl = sapply(symbol2entrez,function(x,y)unique(y[x]),y=entrez2ensmbl)
get_fuma_gtex_deg_data<-function(path,thr=0.1){
  res = read.delim(paste(path,"DEGgeneral.txt",sep=''),stringsAsFactors = F)
  res = res[res$adjP<=thr,]
  gene_sets = lapply(res$genes,function(x)strsplit(x,split=':')[[1]])
  gene_sets = lapply(gene_sets,function(x,y)sapply(x,function(a,b)unlist(b[[a]]),b=ensembl2symbol))
  names(gene_sets) = paste(res$GeneSet,res$Category,sep=';')
  return(list(enrichment=res,gene_sets=gene_sets))
}
rep_analysis_flow<-function(data_2_3,data_2_3_leadsnps,data_1_3){
  # Analysis 1: Lead SNPs
  ind_sig_snps_2_3 = paste(data_2_3_leadsnps$chr,data_2_3_leadsnps$pos,sep=';')
  ind_sig_repdata = data_1_3[ind_sig_snps_2_3,]
  ind_snp_direct_rep = ind_sig_repdata$X13_p<0.01
  
  # Analysis 2: by regions
  ind_snps = data_2_3_leadsnps$rsID
  region2pvals = list()
  for(snp in ind_snps){
    curr_rows = data_2_3[data_2_3$IndSigSNP==snp,]
    curr_snp_pos = paste(curr_rows$chr,curr_rows$pos,sep=';')
    curr_rep_data = data_1_3[curr_snp_pos,]
    ps = curr_rep_data$X13_p
    names(ps) = curr_snp_pos
    region2pvals[[snp]] = ps[!is.na(ps)]
  }
  region2qvals = lapply(region2pvals,p.adjust,method="fdr")
  ind_snp_block_rep = sapply(region2qvals,function(x)any(x<0.1))
  table(ind_snp_block_rep,ind_snp_direct_rep)
  ind_snp_mafs = data_2_3[ind_snps,]$MAF
  m = cbind(ind_sig_snps_2_3,ind_snp_direct_rep,ind_snp_block_rep,ind_snp_mafs,ind_sig_repdata$X13_p)
  colnames(m) = c("position_id","direct_rep","block_rep","MAF","rep_p")
  rownames(m) = ind_snps
  return(m)
}
##############################################
##############################################
##############################################

# Run replication analysis, get the display items
replication_analysis_results = list()
for (nn in names(replication_analysis_files)){
  data_2_3 = get_fuma_all_snp_data(replication_analysis_files[[nn]]$`67`)
  data_2_3_leadsnps = get_fuma_lead_snp_data(replication_analysis_files[[nn]]$`67`)
  data_1_3 = read.delim(replication_analysis_files[[nn]]$rep,
                        row.names = NULL,header=T,stringsAsFactors=F,sep=" ")
  rownames(data_1_3) = paste(data_1_3[,1],data_1_3[,2],sep=';')
  replication_analysis_results[[nn]] = rep_analysis_flow(data_2_3,data_2_3_leadsnps,data_1_3)  
}
rep_numbers = sapply(replication_analysis_results,function(x)apply(x[,2:3],2,function(y)sum(y==TRUE)/length(y)))
rownames(rep_numbers) = c("Direct","Block")
# Figures
barplot(rep_numbers,legend=T,beside=T,args.legend = list(x="topleft"),ylab = "Percent replicated",ylim=c(0,1.1))
exercise_rep_maps = data.frame(Replicated=replication_analysis_results$exercise[,3],MAF=as.numeric(replication_analysis_results$exercise[,4]))
p <- ggplot(exercise_rep_maps, aes(x=Replicated, y=MAF,fill=Replicated)) + geom_violin(scale="width")
p + theme(legend.position="none")

# GWAScatalog summaries
# Direct SNP annotation
default_fuma_leadsnps = get_fuma_lead_snp_data(fuma_pooled_analysis_maffilter_path$exercise)
default_fuma_gwascat = get_fuma_gwascatalog_results(fuma_pooled_analysis_maffilter_path$exercise)
trait2leadsnps = tapply(default_fuma_gwascat$GenomicLocus,default_fuma_gwascat$Trait,unique)
trait2leadsnps[["NA"]] = setdiff(default_fuma_leadsnps$GenomicLocus,default_fuma_gwascat$GenomicLocus)
trait2numsnps = sort(sapply(trait2leadsnps,length),decreasing=T)
names(trait2numsnps) = gsub(names(trait2numsnps) ,pattern="Waist circumference",replace = "Waist circ")
names(trait2numsnps) = gsub(names(trait2numsnps) ,pattern="adjusted",replace = "adj")
names(trait2numsnps) = gsub(names(trait2numsnps) ,pattern="joint analysis ",replace = "")
exercise_trait2numsnps = trait2numsnps
default_fuma_leadsnps = get_fuma_lead_snp_data(fuma_pooled_analysis_maffilter_path$recovery)
default_fuma_gwascat = get_fuma_gwascatalog_results(fuma_pooled_analysis_maffilter_path$recovery)
trait2leadsnps = tapply(default_fuma_gwascat$GenomicLocus,default_fuma_gwascat$Trait,unique)
trait2leadsnps[["NA"]] = setdiff(default_fuma_leadsnps$GenomicLocus,default_fuma_gwascat$GenomicLocus)
trait2numsnps = sort(sapply(trait2leadsnps,length),decreasing=T)
names(trait2numsnps) = gsub(names(trait2numsnps) ,pattern="Waist circumference",replace = "Waist circ")
names(trait2numsnps) = gsub(names(trait2numsnps) ,pattern="adjusted",replace = "adj")
recovery_trait2numsnps = trait2numsnps
recovery_trait2numsnps = recovery_trait2numsnps[recovery_trait2numsnps>1 | names(recovery_trait2numsnps)=="NA"]
exercise_trait2numsnps = exercise_trait2numsnps[exercise_trait2numsnps>1]
par(mar=c(5.1, 20.1, 4.1, 2.1))
rec_cols = rep("cyan",length(recovery_trait2numsnps))
rec_cols[names(recovery_trait2numsnps)=="NA"] = "black"
barplot(recovery_trait2numsnps[1:length(recovery_trait2numsnps)],
        horiz = T,las=2,xlab = "Number of regions",col=rec_cols,density = 80,space=0)
ex_cols = rep("cyan",length(exercise_trait2numsnps))
ex_cols[names(exercise_trait2numsnps)=="NA"] = "black"
barplot(exercise_trait2numsnps[1:length(exercise_trait2numsnps)],
        horiz = T,las=2,xlab = "Number of regions",cex.names=0.8,
        col=ex_cols,density = 80,space=0)

# Useful plots for showing enrichments
pairwise_barplot<-function(x1,x2,ns=NULL,col1="red",col2="blue",cex.text=1.5,cex.axis=1.5,
                           xlab1 = "Proportion of overlapping\ngenes in gene sets",
                           xlab2 = "Log10 adj\nP-value",...){
  if(!is.null(ns)){names(x2) = ns}
  if(all(x2>0)){x2 = -x2}
  ord = order(x2,decreasing = T)
  par(mfrow=c(1,3))
  dummyvals = rep(0,length(x1))
  names(dummyvals) = names(x2[ord])
  par(mar = c(5.1,0,4.1,0))
  xx = barplot(dummyvals,plot=F,beside = T)
  barplot(xx[,1],xlim=c(0, 1), ann = F, bty = 'n',xaxt = 'n', yaxt = 'n',beside=T,horiz=T,col="white",border=F)
  for(i in 1:length(x1)){
    text(x = 0.01, y = xx[i,1], names(x2[ord])[i], 
         cex = cex.text, col = "black", family="serif", font=2, adj=0)
  }
  par(mar=c(5.1, 0, 4.1, 0))
  barplot(x2[ord],horiz = T,space=0,col=col2,main=xlab2,las=2,names.arg = "",
          xlim=c(min(x2),-0.01),xpd = F,cex.axis = cex.axis)
  par(mar = c(5.1,0,4.1,2))
  barplot(x1[ord],horiz = T,space=0,col=col1,main=xlab1,las=2,cex.axis = cex.axis)
}
format_gene_gwascat<-function(x){
  rownames(x) = x$GeneSet
  x = x[,c(3:6)]
  x$Proportion = x[,2]/x[,1]
  x = x[,5:4]
  x[,2] = log(x[,2],10)
  return(x)
}
split_name_by_num_words<-function(s,n){
  arr = strsplit(s,split=" ")[[1]]
  if(n>length(arr)){return(s)}
  news = paste(paste(arr[1:n],collapse=" "),paste(arr[(n+1):length(arr)],collapse=" "),sep="\n")
  return(news)
}

# Gene-based enrichments
gene_enrichments = lapply(fuma_pooled_analysis_maffilter_path,get_fuma_gene_sets_results)
gene_gwascat_enrichments = lapply(gene_enrichments,function(x)x[x[,1]=="GWAScatalog",])
gene_gwascat_enrichments = lapply(gene_gwascat_enrichments,format_gene_gwascat)
# Recovery plot
rownames(gene_gwascat_enrichments[[2]])[2] = 
  split_name_by_num_words(rownames(gene_gwascat_enrichments[[2]])[2],3)
rownames(gene_gwascat_enrichments[[2]])[4] = 
  split_name_by_num_words(rownames(gene_gwascat_enrichments[[2]])[4],3)
pairwise_barplot(gene_gwascat_enrichments[[2]][,1],
                 gene_gwascat_enrichments[[2]][,2],
                 rownames(gene_gwascat_enrichments[[2]]),
                 cex.text=1.5)
# Exercise
rownames(gene_gwascat_enrichments[[1]])[2] = 
  split_name_by_num_words(rownames(gene_gwascat_enrichments[[1]])[2],4)
rownames(gene_gwascat_enrichments[[1]])[9] = 
  split_name_by_num_words(rownames(gene_gwascat_enrichments[[1]])[9],2)
pairwise_barplot(gene_gwascat_enrichments[[1]][1:15,1],
                 gene_gwascat_enrichments[[1]][1:15,2],
                 rownames(gene_gwascat_enrichments[[1]])[1:15],
                 cex.text=1.3)

# The third analysis involves preparation of data for network-based
#   analysis to seek pathways e.g., using ResponseNet
get_pairwise_overlap_and_pval<-function(s1,s2,bg_min_size){
  q = length(intersect(s1,s2))
  m = length(s1)
  n = bg_min_size - m
  k = length(s2)
  pval = phyper(q-1,m,n,k,lower.tail = F)
  return(pval)
}
# Exercise
ex_fuma_gene_data = get_fuma_mapped_gene_data(fuma_pooled_analysis_maffilter_path$exercise)
ex_fuma_gtex_data = get_fuma_gtex_deg_data(fuma_pooled_analysis_maffilter_path$exercise)
# Recovery
rec_fuma_gene_data = get_fuma_mapped_gene_data(fuma_pooled_analysis_maffilter_path$recovery)
rec_fuma_gtex_data = get_fuma_gtex_deg_data(fuma_pooled_analysis_maffilter_path$recovery)
# RHR
rhr_fuma_gene_data = get_fuma_mapped_gene_data(fuma_pooled_analysis_maffilter_path$rhr)
# Meta-analysis data
meta_analysis_gene_data = get(load("../MoTrPAC/PA_database/Exercise_data_analysis_gene_sets.RData"))
meta_analysis_gene_data = meta_analysis_gene_data[grepl(names(meta_analysis_gene_data),pattern="META_")]
meta_analysis_raw_results = get(load("../MoTrPAC/PA_database/tp_meta_analysis_results.RData"))

all_mut_genes = union(ex_fuma_gene_data$symbol,rec_fuma_gene_data$symbol)
all_mut_genes_ids = as.character(union(ex_fuma_gene_data$entrezID,rec_fuma_gene_data$entrezID))
inter_genes = intersect(names(meta_analysis_gene_data$`META_longterm,muscle`),all_mut_genes_ids)
inter_genes_names = unlist(entrez2symbol[inter_genes])
inter_genes_ensembl = entrez2ensmbl[as.character(inter_genes)]

fuma_gene_set1 = ex_fuma_gene_data$symbol
ex_ps = sapply(meta_analysis_gene_data,get_pairwise_overlap_and_pval,s2=fuma_gene_set1,bg_min_size=20000)
fuma_gene_set1 = rec_fuma_gene_data$symbol
rec_ps = sapply(meta_analysis_gene_data,get_pairwise_overlap_and_pval,s2=fuma_gene_set1,bg_min_size=20000)
fuma_gene_set1 = rhr_fuma_gene_data$symbol
rhr_ps = sapply(meta_analysis_gene_data,get_pairwise_overlap_and_pval,s2=fuma_gene_set1,bg_min_size=20000)
meta_analysis_overlap_ps = c(ex_ps[3],rec_ps[3],rhr_ps[3])
names(meta_analysis_overlap_ps) = c("ExerciseHR","Recovery","RHR")
par(mfrow=c(1,1),mar=c(4,4,4,4))
barplot(-log(meta_analysis_overlap_ps,10),main="Significance of overlap with meta-analysis",
        ylab = "-log10 P-value",space = 0)

# Summarize the results in one big table
m = cbind(is.element(all_mut_genes,set=ex_fuma_gene_data$symbol),is.element(all_mut_genes,set=rec_fuma_gene_data$symbol))
rownames(m) = all_mut_genes
colnames(m) = c("ExerciseHR","Recovery")
all_expr_data = meta_analysis_raw_results$`longterm,muscle`
rownames(all_expr_data) = entrez2symbol[rownames(all_expr_data)]
all_genes = union(all_mut_genes,meta_analysis_gene_data$`META_longterm,muscle`)
all_expr_data = all_expr_data[meta_analysis_gene_data$`META_longterm,muscle`,]
all_expr_data = all_expr_data[,1:3]
mut_deg_matrix = matrix(NA,ncol=5,nrow=length(all_genes))
rownames(mut_deg_matrix) = all_genes
mut_deg_matrix[rownames(all_expr_data),1:3] = all_expr_data
colnames(mut_deg_matrix) = c(colnames(all_expr_data),"MutExerciseHR","MutRecovery")
mut_deg_matrix[rownames(m),4:5] = m
mut_deg_matrix = cbind(mut_deg_matrix,is.element(rownames(mut_deg_matrix),set=ex_fuma_gtex_data[[2]]$`Muscle;DEG.up`))
colnames(mut_deg_matrix)[6] = "GTeXMuscle"
gene_type = rep("DEG",nrow(mut_deg_matrix))
gene_type[mut_deg_matrix[,4]==1 & mut_deg_matrix[,5]==1] = "MutBoth"
gene_type[mut_deg_matrix[,4]==1 & mut_deg_matrix[,5]!=1] = "MutExercise"
gene_type[mut_deg_matrix[,4]!=1 & mut_deg_matrix[,5]==1] = "MutRecovery"
mut_deg_matrix = cbind(mut_deg_matrix,gene_type)
mut_deg_matrix = mut_deg_matrix[,-c(4:5)]
mut_deg_matrix[inter_genes_names,"gene_type"] = "MutatedAndDEG"
write.table(mut_deg_matrix,file="mutation_deg_data_for_network_analysis.txt",sep="\t",quote=F)

# # Create input files for ResponseNet
# # Targets: DEGs
# tars = t(t(meta_analysis_gene_data$`META_longterm,muscle`))
# tars = unlist(entrez2ensmbl[rownames(tars)])
# tars = tars[!is.na(tars)]
# tars = t(t(tars));rownames(tars) = tars[,1]
# tars = cbind(rownames(tars),rep("Gene",length(tars)))
# write.table(file="gwas_interpretation/responsenet/pooled_maffilter_targets.txt",tars,row.names = F,quote=F,col.names = F,sep="\t")
# 
# # sources: mutated genes, target: expression genes from meta-analysis
# sources = t(t(all_mut_genes_ids))
# sources = unlist(entrez2ensmbl[sources[,1]])
# sources = sources[!is.na(sources)]
# sources = t(t(sources))
# s_weights = rep(0.5,length(sources))
# s_weights[is.element(rownames(sources),set=rownames(tars))] = 1
# sources = cbind(sources,rep("Gene",length(sources)),s_weights)
# write.table(file="gwas_interpretation/responsenet/pooled_maffilter_sources.txt",sources,row.names = F,quote=F,col.names = F,sep="\t")
# 
# # Create input files for ResponseNet
# # sources: ONLY exercise mutated genes, target: expression genes from meta-analysis
# sources = t(t(as.character(ex_fuma_gene_data$entrezID)))
# sources = unlist(entrez2ensmbl[sources[,1]])
# sources = sources[!is.na(sources)]
# sources = t(t(sources))
# s_weights = rep(0.5,length(sources))
# s_weights[is.element(rownames(sources),set=rownames(tars))] = 1
# sources = cbind(sources,rep("Gene",length(sources)),s_weights)
# write.table(file="gwas_interpretation/responsenet/pooled_maffilter_exercise_only_sources.txt",sources,row.names = F,quote=F,col.names = F,sep="\t")

# Methods for comparing a pair of GWAS runs in different ways
library(Vennerable)
region_based_comparison<-function(path1,path2){
  res1 = get_fuma_lead_snp_data(path1)
  res2 = get_fuma_lead_snp_data(path2)
  res1_allsnps = unlist(sapply(res1$IndSigSNPs,function(x)strsplit(x,split=';')[[1]]))
  res2_allsnps = unlist(sapply(res2$IndSigSNPs,function(x)strsplit(x,split=';')[[1]]))
  V = Venn(list(res1_allsnps,res2_allsnps))
  plot(V,doWeights=T)
}
snp_based_comparison<-function(path1,path2,ukbb_maf = snp2maf){
  res1 = get_fuma_all_snp_data(path1)
  res1 = res1[!is.na(res1$gwasP),]
  res2 = get_fuma_all_snp_data(path2)
  res2 = res2[!is.na(res2$gwasP),]
  all_snps = union(res1$rsID,res2$rsID)
  summat = cbind(is.element(all_snps,set=res1$rsID),is.element(all_snps,set=res2$rsID))
  fuma_mafs = rep(0,length(all_snps));names(fuma_mafs) = all_snps
  fuma_mafs[res1$rsID] = res1$MAF
  fuma_mafs[res2$rsID] = res2$MAF
  ukbb_mafs = ukbb_maf[all_snps]
  plot(fuma_mafs,ukbb_mafs);abline(0,1)
  rownames(summat) = all_snps
  summat = cbind(summat,ukbb_mafs,fuma_mafs)
  colnames(summat) = c("in_set1","in_set2","MAF_UKBB","MAF_1000G")
  return(summat)
}
enrichment_based_comparison<-function(path1,path2,types=c("GWAScatalog"),adjP_thr=0.01){
  res1 = get_fuma_gene_sets_results(path1)
  res2 = get_fuma_gene_sets_results(path2)
  res1 = res1[is.element(res1$Category,set=types),]
  res2 = res2[is.element(res2$Category,set=types),]
  res1 = res1[res1$adjP<=adjP_thr,]
  res2 = res2[res2$adjP<=adjP_thr,]
  V = Venn(list(s1 = res1$GeneSet,s2=res2$GeneSet))
  plot(V,doWeights=T)
}

# Saturation analysis
num_loci = sapply(saturation_analysis_paths$exercise,function(x)length(unique(get_fuma_lead_snp_data(x)$GenomicLocus)))
num_leads = sapply(saturation_analysis_paths$exercise,function(x)length(unique(get_fuma_lead_snp_data(x)$rsID)))
par(mfrow=c(1,1),mar=c(5,5,5,5))
plot(x=names(num_leads),y=num_loci,type="l",pch=20,lwd=4,ylab = "Number of regions",col="darkblue",
     xlab = "Percent sampled, Exercise HR",ylim=c(0,250),xlim = c(30,105))
plot(x=names(num_leads),y=num_loci,type="l",pch=20,lwd=4,ylab = "Number of lead SNPs",col="darkblue",
     xlab = "Percent sampled, Exercise HR",xlim = c(30,105))

# Exercise HR: Get pairwise comparison between 67% and 85%
p1 = saturation_analysis_paths$exercise$`67`
p2 = saturation_analysis_paths$exercise$`85`
snp_comp = snp_based_comparison(p1,p2,snp2maf)
snp_sets = rep("85",nrow(snp_comp))
snp_sets[snp_comp[,1]==1] = "67"
snp_sets[snp_comp[,1]==1 & snp_comp[,2]==1] = "both"

s1 = rownames(snp_comp)[snp_comp[,1]==1]
s2 = rownames(snp_comp)[snp_comp[,2]==1]
V = Venn(list('67'=s1,'85'=s2))
gp = VennThemes(compute.Venn(V))
gp_cols = sapply(gp$Face,function(x)x$col)
plot(V,doWeights=T)

vioplot_data = data.frame(list(MAF=snp_comp[,3],set=snp_sets))
vioplot_data[,1] = as.numeric(as.character(vioplot_data[,1]))
vioplot_data[,2] = factor(vioplot_data[,2],levels = c("67","both","85"))
p <- ggplot(vioplot_data, aes(x=set, y=MAF,fill=set)) + geom_violin(scale="width")
p + scale_fill_manual(values=c('67'= unname(gp_cols['10.10']), '85'=unname(gp_cols['01.01']), "both"=unname(gp_cols[1]))) + theme_classic()

# Exercise HR: Stratify p-values by the MAF
# raw_data = read.delim("gwas_interpretation/fuma_in_files/HR_fitnes_raw.txt",sep=" ",stringsAsFactors = F,header = T)
# snp_run_ids = paste(raw_data[,1],raw_data[,2],sep=';')
# snp_id_data = read.delim("gwas_interpretation/id.txt",sep=" ",stringsAsFactors = F,header = T)
# snp_run_ids2 = paste(snp_id_data[,1],snp_id_data[,2],sep=';')
# snp_run_ids2_table = table(snp_run_ids2)
# snps_with_dups = names(which(snp_run_ids2_table>1))
# dup_inds = is.element(snp_run_ids2,snps_with_dups)
# loc2rsid = snp_id_data[!dup_inds,3]
# names(loc2rsid) = snp_run_ids2[!dup_inds]
# dup_data =  snp_id_data[dup_inds,]
# raw_data = cbind(raw_data,loc2rsid[snp_run_ids])
# curr_mafs = snp2maf[as.character(raw_data[,4])]
# high_maf_rows = curr_mafs > 0.005
# table(high_maf_rows)
# ps1 = raw_data[high_maf_rows,3]
# ps2 = raw_data[!high_maf_rows,3]
# samp1 = sample(ps1)[1:100000]
# samp2 = sample(ps2)[1:100000]
# par(mfrow=c(2,1))
# qqplot(-log(runif(10000),10),-log(samp1,10),main = "MAF > 0.5%",
#        xlab="Expected -log10 P-value",ylab = "Observed -log10 P-value",cex=0.8);abline(0,1)
# qqplot(-log(runif(10000),10),-log(samp2,10),main="MAF <= 0.5%",
#        xlab="Expected -log10 P-value",ylab = "Observed -log10 P-value",cex=0.8);abline(0,1)
# qqplot(samp1,samp2);abline(0,1)

# # Look at the hemoglobin SNPs
# # These selected snps had p<1e-4
# # This computation was done directly using the MR code
# hemo_lead_snps = c("rs9705925","rs7299011","rs13321817","rs12648115","rs9784505","rs16891982","rs3756298","rs1896910" )
# exercise_leadsnps = get_fuma_lead_snp_data(fuma_pooled_analysis_maffilter_path$exercise)
# recovery_leadsnps = get_fuma_lead_snp_data(fuma_pooled_analysis_maffilter_path$recovery)
# table(is.element(hemo_lead_snps,set=exercise_leadsnps$rsID))
# table(is.element(hemo_lead_snps,set=recovery_leadsnps$rsID))
# exercise_snps = get_fuma_all_snp_data(fuma_pooled_analysis_maffilter_path$exercise)
# write.table(exercise_snps[hemo_lead_snps,c("MAF","nearestGene")],sep="\t",quote=F)

