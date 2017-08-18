# On Sherlock
geno_path = '/oak/stanford/groups/euan/projects/ukbb/gwas/snp_extract/fitness'
pheno_path = ''
icd_path = '/oak/stanford/groups/euan/projects/ukbb/code/anna_code/icd/icd_matrix.txt'

# Local
geno_path = '/Users/David/Desktop/ukbb/mr/genotypes/fitness/'
pheno_path = '/Users/David/Desktop/ukbb/fitness_analysis_final_fitness_scores.RData'
icd_path = '/Users/David/Desktop/ukbb/mr/icd_matrix.txt'
snps_fitness_table = '/Users/David/Desktop/ukbb/top.snps.xlsx'

# read the exposure data
get(load(pheno_path))
expo = fitness_scores_matrix

# read the genotypes
geno_files = list.files(geno_path)
geno_files = geno_files[grepl("\\.raw$",geno_files)]
geno_data = NULL
for(f in geno_files){
  print(f)
  snp_data = read.delim(paste(geno_path,f,sep=''),header=T,sep=" ")
  rownames(snp_data) = snp_data[,1]
  snp_data = snp_data[,-c(1:6)]
  snp_data = snp_data[,!grepl("_HET$",colnames(snp_data))]
  if(!is.null(geno_data) && any(rownames(geno_data)!=rownames(snp_data))){
    print ("ERROR: row names do not match")
    break
  }
  if(is.null(geno_data)){
    geno_data = snp_data
    next
  }
  geno_data = cbind(geno_data,snp_data)
}
colnames(geno_data) = gsub(colnames(geno_data),pattern = "_.$",replace="")
save(geno_data,file=paste(geno_path,"fitness_geno_data_from_raw_files.RData",sep=''))

# read the outcome
icd_data = read.delim(icd_path)
rownames(icd_data) = icd_data[,1]
heart_disease_codes = c(
  paste("I",21:22,sep=""), #myocardial infraction
  "I42","I50", #heart failure
  "I48", #atrial fibrilation
  "I471" #supraventricular tachycardia 
)
hypertension_codes = paste("I",10:15,sep="")
diabetes_codes = paste("IE",10:14,sep="")

heart_disease_outcome = as.numeric(apply(icd_data[,intersect(colnames(icd_data),heart_disease_codes)],1,any))
diabetes_outcome = as.numeric(icd_data[,intersect(colnames(icd_data),hypertension_codes)])
hypertension_outcome = as.numeric(apply(icd_data[,intersect(colnames(icd_data),diabetes_codes)],1,any))
names(heart_disease_outcome) = rownames(icd_data)
names(diabetes_outcome) = rownames(icd_data)
names(hypertension_outcome) = rownames(icd_data)

# Read the snp lists
library(xlsx)
fitness_snp_data = as.matrix(read.xlsx2(snps_fitness_table,sheetIndex = 1))
fitness_snps = fitness_snp_data[fitness_snp_data[,"ancestry"]=="european" &
    (grepl(fitness_snp_data[,"trait"],pattern="Predicted_HR") | grepl(fitness_snp_data[,"trait"],pattern="Rest"))& 
    fitness_snp_data[,"score_adjustment"] == "conservative"
    ,"snp"]
fitness_snps = fitness_snp_data[fitness_snp_data[,"ancestry"]=="european" &
    (grepl(fitness_snp_data[,"trait"],pattern="Predicted_HR") | grepl(fitness_snp_data[,"trait"],pattern="Rest")) & 
    fitness_snp_data[,"score_adjustment"] == "simple"
    ,"snp"]

# Format all data to be on the same subjects
subjs = rownames(expo)[!is.na(expo[,2])]
subjs = intersect(subjs,rownames(geno_data))
expo_v = expo[subjs,2]
outcome_v = heart_disease_outcome[subjs]
snps = geno_data[subjs,intersect(fitness_snps,colnames(geno_data))]
snps = as.matrix(snps)
dim(snps)
gc()

# get betas and their sds, also merge snps based on betas
lm_expo_vs_g = lm(expo_v~snps)
bx_obj = summary(lm_expo_vs_g)$coefficients
bx = bx_obj[-1,1]
bxse = bx_obj[-1,2]
lm_outcome_vs_g = lm(outcome_v~snps)
by_obj = summary(lm_outcome_vs_g)$coefficients
by = by_obj[-1,1]
byse = by_obj[-1,2]

# merge snps by beta and rerun
weighted_causal_v = snps %*% bx
lm_expo_vs_g = lm(expo_v~weighted_causal_v)
bx_obj = summary(lm_expo_vs_g)$coefficients
bx = bx_obj[-1,1]
bxse = bx_obj[-1,2]
lm_outcome_vs_g = lm(outcome_v~weighted_causal_v)
by_obj = summary(lm_outcome_vs_g)$coefficients
by = by_obj[-1,1]
byse = by_obj[-1,2]

# Perform MR
library('MendelianRandomization')
mr_in = mr_input(bx,bxse,by,byse)
mr_plot(mr_in)
mr_maxlik(mr_in)
mr_ivw(mr_in)
# multivariate analysis
mr_median(mr_in)
mr_egger(mr_in)
mr_allmethods(mr_in)


# Get the data from the ped and map files (instead of raw)
library(HardyWeinberg)
transform_geno<-function(v,return_numeric=T){
  v = as.matrix(v)
  major = names(sort(table(v[,1]),decreasing=T))[1]
  v[v==major]="A"
  v[v!="A" & v!="0"]="B"
  v[v=="0"] = "O"
  ab = apply(v,1,paste,collapse='')
  ab[ab=="BA"] = "AB"
  #hw_test = HWLratio(table(ab)[c("AA","AB","BB")])
  if(return_numeric){
    ab[ab=="AA"] = "0"
    ab[ab == "AB"] = "1"
    ab[ab=="BB"] = "2"
    ab[ab=="OO"] = NA
    return(as.numeric(ab))
  }
  return(ab)
}

# read the genotype data, get the current subjects at a time and merge the genotypes
library(Matrix)
geno_files = list.files(geno_path)
geno_files = geno_files[!grepl(geno_files,pattern="RData$")]
geno_names = unique(sapply(geno_files,function(x)paste(strsplit(x,split='\\.')[[1]][1:2],collapse='.')))
geno_parsed_data = NULL
for(g in geno_names){
  print(g)
  snp_data = read.delim(paste(geno_path,g,'.map',sep=''),header=F)
  rdata = read.delim(paste(geno_path,g,'.ped',sep=''),sep=" ",header = F)
  colnames(rdata) = c("famid", "pid", "fatid", "motid", "sex","affected",sapply(as.character(snp_data[,2]),rep,times=2))
  if(ncol(rdata) != 6 + 2*nrow(snp_data)){
    print ("ERROR: number of snps does not match the ped data")
    break
  }
  gdata = c()
  for (snp in as.character(snp_data[,2])){
    inds = colnames(rdata)==snp
    gv = transform_geno(rdata[,inds])
    gdata = cbind(gdata,gv)
    colnames(gdata)[ncol(gdata)] = snp
  }
  rownames(gdata) = rdata[,1]
  gdata = Matrix(gdata)
  if(!is.null(geno_parsed_data) && any(rownames(geno_parsed_data)!=rownames(gdata))){
    print("ERROR:rownames do not match")
    break
  }
  if(is.null(geno_parsed_data)){
    geno_parsed_data = gdata
  }
  else{
    geno_parsed_data = cbind(geno_parsed_data,gdata)
  }
  gc()
}
geno_data = geno_parsed_data
save(geno_data,file=paste(geno_path,"fitness_geno_data_from_ped_files.RData",sep=''))
