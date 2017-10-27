# Functional analysis of the snp results file
try({setwd("/Users/david/Desktop/ukbb/")})
library(xlsx)

###############################################
###############################################
############# Helper functions ################
###############################################
###############################################

extract_genes_GREAT<-function(x){
  arr = strsplit(x,split = " ")[[1]]
  arr = arr[!grepl('\\(',arr)]
  return(arr)
}
get_pairwise_matrix<-function(l,func){
  n = length(l)
  m = matrix(NA,nrow=n,ncol=n)
  rownames(m) = names(l)
  colnames(m) = names(l)
  for(i in 2:n){
    for(j in 1:(i-1)){
      m[i,j] = func(l[[i]],l[[j]])
      m[j,i]=m[i,j]
    }
  }
  return(m)
}
get_unique_in_ind<-function(j,l){
  ll = l[-j]
  return(setdiff(l[[j]],unique(unlist(ll))))
}
# xx is a matrix with the columns: "snp", "chr", "pos"
print_GREAT_input<-function(xx,fname){
  all_snps_bed = unique(xx[,c("snp","chr","pos")])
  all_snps_bed = cbind(all_snps_bed,as.numeric(as.character(all_snps_bed[,"pos"]))+1)
  colnames(all_snps_bed) = c("name","chrom","chromStart","chromEnd")
  all_snps_bed = all_snps_bed[,c(2:4,1)]
  all_snps_bed[,1] = paste("chr",all_snps_bed[,1],sep="")
  write.table(all_snps_bed,file=fname,
              sep="\t",row.names=F,col.names = F,quote=F)
}
print_dummy_FUMA_input<-function(xx,fname){
  xx = xx[,c("chr","pos","snp")]
  xx[,1] = paste("chr",xx[,1],sep="")
  xx = unique(xx)
  fuma_in = cbind(xx[,c("chr","pos")],"0.000000000001")
  fuma_in = cbind(fuma_in,xx[,"snp"])
  colnames(fuma_in) = c("chromosome","position","P-value","rsID")
  write.table(fuma_in,file=fname,sep="\t",row.names=F,col.names = T,quote=F)
  
}
list_merge<-function(l1,l2){
  l = list()
  for(nn in union(names(l1),names(l2))){
    l[[nn]] = union(get_from_list(l1,nn),get_from_list(l2,nn))
  }
  return(l)
}
get_from_list<-function(l,nn,val=NULL){
  if(nn %in% names(l)){
    return (l[[nn]])
  }
  return(val)
}
############################

###############################################
###############################################
#################### End ######################
###############################################
###############################################

###############################################
###############################################
#### Aug 2017: Load the data, print files #####
###############################################
###############################################

dirr = 'gwas/results_august_2017/'
dirr = 'comp/results/simple_norm/'
files = list.files(dirr)
files = files[grepl('\\.txt$',files)]
f_matrices = list()
for (f in files){
  f_matrices[[f]] = read.delim(paste(dirr,f,sep=''),sep=' ',header=F)
  f_matrices[[f]] = f_matrices[[f]][!is.na(f_matrices[[f]][,1]),]
  colnames(f_matrices[[f]]) = c("chr","snp","pos","p")
}

dir.create('gwas/results_august_2017/fuma_input')
for (f in files){
  fuma_in = print_dummy_FUMA_input(f_matrices[[f]],paste("gwas/results_august_2017/fuma_input/",f,sep=''))
}

###############################################
###############################################
#################### End ######################
###############################################
###############################################

###############################################
###############################################
####July 2017: Load the data, print files #####
###############################################
###############################################

raw = as.matrix(read.xlsx2("top.snps.xlsx",sheetIndex = 1))
forlists = raw[,c("snp","score_adjustment","ancestry","chr","pos","freq")]
# Partition the data frame by traits and their correction
# Also, merge the different execise/non exercise traits
trait2data = list()
exercise_regex = c("Predicted","Max_achi","slopes","Rest_HR")
for(tr in unique(raw[,"trait"])){
  trname = gsub(tr,pattern='normalized_ukb_imp',replace='')
  trait2data[[trname]]=forlists[raw[,"trait"]==tr,-2]
  if(any(sapply(exercise_regex,grepl,x=tr)) && grepl("euro",tr)){
    if(grepl("conserv",tr)){
      trait2data[["Exercise_conservative_euro"]] = rbind(trait2data[["Exercise_conservative_euro"]],trait2data[[trname]])
    }
    else{
      trait2data[["Exercise_simple_euro"]] = rbind(trait2data[["Exercise_simple_euro"]],trait2data[[trname]])
    }
  }
  else{
    if(grepl("conserv",tr)){
      trait2data[["Non exercise_conservative_euro"]] = rbind(trait2data[["Non exercise_conservative_euro"]],trait2data[[trname]])
    }
    else{
      trait2data[["Non exercise_simple_euro"]] = rbind(trait2data[["Non exercise_simple_euro"]],trait2data[[trname]])
    }
  }
}
for(nn in names(trait2data)){
  if(is.null(dim(trait2data[[nn]]))){
    currnames = names(trait2data[[nn]])
    trait2data[[nn]] = matrix(trait2data[[nn]],nrow=1)
    colnames(trait2data[[nn]]) = currnames
  }
  else{
    print(dim(trait2data[[nn]]))
    trait2data[[nn]] = unique(trait2data[[nn]])
    print(dim(trait2data[[nn]]))
  }
}
sort(sapply(trait2data,nrow))

# Create input files for the online tools:
# a BED file for GREAT and another file for FUMA
all_snps_bed = unique(raw[,c("snp","chr","pos")])
all_snps_bed = cbind(all_snps_bed,as.numeric(as.character(all_snps_bed[,"pos"]))+1)
colnames(all_snps_bed) = c("name","chrom","chromStart","chromEnd")
all_snps_bed = all_snps_bed[,c(2:4,1)]
all_snps_bed[,1] = paste("chr",all_snps_bed[,1],sep="")
# Create input files for FUMA and GREAT
great_dir = "GREAT_input_files2/"
fuma_dir = "FUMA_input_files2/"
dir.create(great_dir);dir.create(fuma_dir)
# For GREAT:
print_GREAT_input(raw,fname=paste(great_dir,"all_snps_GREAT.bed",sep=''))
# For FUMA:
print_dummy_FUMA_input(raw,fname=paste(fuma_dir,"all_snps_FUMA.bed",sep=''))
for(nn in names(trait2data)){
  mat = trait2data[[nn]]
  nn = gsub(nn,pattern = " ",replace="_")
  nn = gsub(nn,pattern = "=",replace="_")
  if(nrow(mat)>1){
    print_GREAT_input(mat,fname=paste(great_dir,nn,".bed",sep=''))
    print_dummy_FUMA_input(mat,fname=paste(fuma_dir,nn,".bed",sep=''))
  }
}

###############################################
###############################################
#################### End ######################
###############################################
###############################################

###############################################
###############################################
########### Load snp2gene, analyze ############
###############################################
###############################################

# Analyze using GREAT results
# load snp2gene file: a mapping of snp ids to a list of genes
snp2gene_raw = as.matrix(read.delim("GREAT_snp2gene.txt"))[,1:2]
snp2gene_lists = sapply(snp2gene_raw[,2],extract_genes_GREAT)
names(snp2gene_lists) = snp2gene_raw[,1]
great_snp2gene_lists = snp2gene_lists

# Analyze using FUMA results
# load snp2gene file: a mapping of snp ids to a list of genes
# genes.txt
fuma_snp2gene_raw = as.matrix(read.delim("FUMA_snp2gene.txt"))
fuma_snp2gene_lists = list()
for (i in 1:nrow(fuma_snp2gene_raw)){
  arr = fuma_snp2gene_raw[i,]
  g = arr["symbol"]
  snps = arr["IndSigSNPs"]
  snps = strsplit(snps,split=':')[[1]]
  for (snp in snps){
    fuma_snp2gene_lists[[snp]] = union(fuma_snp2gene_lists[[snp]],g)
  }
}
# add data from snps.txt
fuma_snp2gene_raw2 = as.matrix(read.delim("FUMA_snp2gene2.txt"))
for (i in 1:nrow(fuma_snp2gene_raw2)){
  arr = fuma_snp2gene_raw2[i,]
  gs = arr["nearestGene"]
  snp = arr["rsID"]
  gs = strsplit(gs,split=':')[[1]]
  for (g in gs){
    fuma_snp2gene_lists[[snp]] = union(fuma_snp2gene_lists[[snp]],g)
  }
}

trait2genes1 = sapply(trait2data,function(x,y)unlist(y[x[,1]]),y=great_snp2gene_lists)
trait2genes2 = sapply(trait2data,function(x,y)unlist(y[x[,1]]),y=fuma_snp2gene_lists)
trait2genes_union = list()
for(nn in union(names(trait2genes1),names(trait2genes2))){
  trait2genes_union[[nn]] = union(trait2genes1[[nn]],trait2genes2[[nn]])
}

###############################################
###############################################
#################### End ######################
###############################################
###############################################

###############################################
###############################################
############## Display items ##################
###############################################
###############################################

# Compare GREAT and FUMA
# shared_snps = intersect(names(great_snp2gene_lists),names(fuma_snp2gene_lists))
# length(shared_snps)
# snp_jaccards = c()
# for(snp in shared_snps){
#   int_size = length(intersect(fuma_snp2gene_lists[[snp]],great_snp2gene_lists[[snp]]))
#   u_size = length(union(fuma_snp2gene_lists[[snp]],great_snp2gene_lists[[snp]]))
#   snp_jaccards[[snp]] = int_size/u_size
# }
# hist(snp_jaccards)
# # get stats
# length(unique(raw$snp))
# length(great_snp2gene_lists)
# length(fuma_snp2gene_lists)
# length(shared_snps)

# trait to total num samples
tt = get(load("physical_fitness_scores_for_GWAS.RData"))
trait2sample_size = apply(tt,2,function(x)sum(!is.na(x)))
tt = get(load("additional_scores_for_GWAS.RData"))
trait2sample_size = c(trait2sample_size,apply(tt,2,function(x)sum(!is.na(x))))
names(trait2sample_size) = gsub(names(trait2sample_size),pattern=" ",replace="_")

# Plots and stats for the presentation
mafs = as.numeric(as.character(raw[,"freq"]))
trait2num_snps = sapply(trait2data,nrow)
eu_num_snps = sort(trait2num_snps[grepl(names(trait2num_snps),pattern="euro")])
names(eu_num_snps) = gsub(names(eu_num_snps),pattern = "_euro",replace="")
par(mar=c(5.1, 15, 1.1, 2.1))
barplot(eu_num_snps,las=2,horiz = T,xlab="Number of SNPs")
eu_num_genes = sapply(trait2genes_union,length)
eu_num_genes = sort(eu_num_genes[grepl(names(eu_num_genes),pattern="euro")])
names(eu_num_genes) = gsub(names(eu_num_genes),pattern = "_euro",replace="")
par(mar=c(5.1, 20, 1.1, 2.1))
all_nums = cbind(eu_num_genes,eu_num_snps[names(eu_num_genes)])
colnames(all_nums) = c("Genes","SNPs")
barplot(t(all_nums),las=2,horiz = T,
        xlab="Number of mapped genes",beside = T,legend=T,
        args.legend = list(x="bottomright",cex=1.2))
eu_trait2data = trait2data[grepl(names(trait2data),pattern="euro")]
names(eu_trait2data) = gsub(names(eu_trait2data),pattern = "_euro",replace="")
mafs_by_trait = lapply(eu_trait2data,function(x)as.numeric(as.character(x[,"freq"])))
par(mar=c(20, 5, 1.1, 2.1))
boxplot(mafs_by_trait,las=2)
# pairwise plot with sample size
sample_sizes = trait2sample_size[rownames(all_nums)]
inds = !is.na(sample_sizes)
par(mfrow=c(1,2))
par(mar=c(6,13,2,0))
names(sample_sizes) = gsub(names(sample_sizes),pattern="conservative",replace="cons")
names(sample_sizes) = gsub(names(sample_sizes),pattern="_after",replace="")
names(sample_sizes) = gsub(names(sample_sizes),pattern="_at_",replace="_")
barplot(-sample_sizes[inds],las=2,horiz = T,
        xlab="Number of subjects",beside = T,legend=F,col="blue",cex.names=0.95,axes=F)
ax_at = seq(0,max(sample_sizes,na.rm=T),length.out = 6)
ax_at = ax_at - ax_at%%10000
ax_at_names = paste(ax_at/10000,"e5",sep="")
ax_at_names[1] = ""
axis(1,at=-ax_at,labels=ax_at_names,las=2)
par(mar=c(6,0,2,2))
barplot(t(all_nums[inds,]),las=2,horiz = T,
        xlab="Set size",beside = T,legend=T,
        args.legend = list(x="bottomright",cex=1.2),names.arg=rep("",sum(inds)))

# print snp lists to the pptx
predicted_hr_cons_snps = trait2snps$`Predicted_HR_at_WD=100_conservative_euro`
print_snps_and_genes(predicted_hr_cons_snps,great_snp2gene_lists,fuma_snp2gene_lists,"Exercise_scores_snp_lists.xlsx","predicted_hr_cons")
predicted_hr_simple_snps = trait2snps$`Predicted_HR_at_WD=100_simple_euro`
print_snps_and_genes(predicted_hr_simple_snps,great_snp2gene_lists,fuma_snp2gene_lists,"Exercise_scores_snp_lists.xlsx","predicted_hr_simple")
rest_cons_snps = trait2snps$Rest_HR_ratios_after_60sec_conservative_euro
print_snps_and_genes(rest_cons_snps,great_snp2gene_lists,fuma_snp2gene_lists,"Exercise_scores_snp_lists.xlsx","rest_cons")
rest_simple_snps = trait2snps$Rest_HR_ratios_after_60sec_simple_euro
print_snps_and_genes(rest_simple_snps,great_snp2gene_lists,fuma_snp2gene_lists,"Exercise_scores_snp_lists.xlsx","rest_simple")


print_snps_and_genes<-function(snps,great_snp2gene_lists,fuma_snp2gene_lists,fname,sheet){
  genes1 = great_snp2gene_lists[snps]
  genes2 = fuma_snp2gene_lists[snps]
  genes = list_merge(genes1,genes2)
  genes = sapply(genes,paste,collapse=',')
  write.xlsx(t(t(genes)),file=fname,sheetName = sheet,append=T)
}

# Look at the exercise snps p-value in the pulse rate analysis
exercise_snps_pr_simple_pval = unique(read.delim("our_snps_pulse_rate_pval_simple.txt",header=F,sep=" "))
rownames(exercise_snps_pr_simple_pval) = exercise_snps_pr_simple_pval[,2]
exercise_snps_pr_cons_pval = unique(read.delim("our_snps_pulse_rate_pval_cons.txt",header=F,sep=" "))
rownames(exercise_snps_pr_cons_pval) = exercise_snps_pr_cons_pval[,2]
snps1 = trait2snps$Exercise_simple_euro
snps2 = trait2snps$Exercise_conservative_euro
ps = list()
ps[["simple, simple"]] = exercise_snps_pr_simple_pval[snps1,9]
snps1[is.na(ps[["simple, simple"]])]
ps[["cons, simple"]] = exercise_snps_pr_simple_pval[snps2,9]
ps[["simple, cons"]] = exercise_snps_pr_cons_pval[snps1,9]
ps[["cons, cons"]] = exercise_snps_pr_cons_pval[snps2,9]
ps = lapply(ps,function(x)x[!is.na(x)])
sapply(ps,median)
sapply(ps,function(x)sum(x<5e-8)/length(x))
sapply(ps,function(x)sum(x>0.01)/length(x))
sapply(ps,function(x)which(x>0.001))

# Print our snps for additional analysis by David H and Anna
our_snps = union(trait2snps$Exercise_conservative_euro,trait2snps$Exercise_simple_euro)
m = matrix("No",nrow=length(our_snps),ncol=2)
rownames(m) = our_snps
m[trait2snps$Exercise_conservative_euro,1]="Yes"
m[trait2snps$Exercise_simple_euro,2] = "Yes"
colnames(m) = c("In_conservative","In_simple")
write.table(m,file="fitness_geas_selected_snps_union.txt",sep="\t",quote=F,row.names = T,col.names = T)


