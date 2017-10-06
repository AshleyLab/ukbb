
# 2/3 vs. 1/3 validation for HR Fitness
RES_DIR = "/oak/stanford/groups/euan/projects/ukbb/da_dh/interpretation/sep18_23_250_1_5_8/simpler_norm"
REP_DIR = "~/work/gwas/sep18_gwas_13/simpler_norm/HR_fitnes_"
PVAL_COL_IN_CHECK=11
HAS_HEADER=F

# 2/3 vs. 80% for HR Fitness
RES_DIR = "/oak/stanford/groups/euan/projects/ukbb/da_dh/interpretation/sep18_23_250_1_5_8/simpler_norm"
REP_DIR = "/oak/stanford/groups/euan/projects/ukbb/da_dh/gwas/sep18_gwas_8/simpler_norm/HR_fitnes_"

# 80% vs. 2/3 for HR Fitness
RES_DIR = "/oak/stanford/groups/euan/projects/ukbb/da_dh/interpretation/sep18_23_250_1_5_8/simpler_norm"
REP_DIR = "/oak/stanford/groups/euan/projects/ukbb/da_dh/gwas/sep18_gwas_8/simpler_norm/HR_fitnes_"

# Other constants
PVAL_COL_IN_CHECK=11
HAS_HEADER=T

setwd(RES_DIR)
fnames = list.files()
fnames = fnames[grep("clumped", fnames)]
result_list = list()
for (f in fnames){
  temp = read.table(f, header = T)
  print(paste(f,"has",nrow(temp),"leading snps"))
  # CHANGE: header = F
  check = read.delim(paste0(REP_DIR, as.character(temp[1,1]),".RHR.glm.linear"), header = HAS_HEADER)
  for (i in 1:nrow(temp)){
    print(i)
    value = F
	# CHANGE: use column names instead of numbers, when possible
    snp_name = as.character(temp[i,"SNP"])
	curr_snps = strsplit(split=",",as.character(temp[i,"SP2"]))[[1]]
	curr_snps = sapply(curr_snps,gsub,pattern="\\(1\\)",replace="")
	curr_snps = union(curr_snps,snp_name)
	curr_snps = setdiff(curr_snps,"NONE")
	inds_in_check = is.element(check[,3],set=curr_snps)
	curr_ps = as.numeric(as.character(check[inds_in_check,PVAL_COL_IN_CHECK]))
	print(curr_ps)
	names(curr_ps) = check[inds_in_check,3]
	result_list[[snp_name]] = curr_ps
  }
}

# get replication scores
p_thr = 0.05
sum(sapply(result_list,function(x)any(p.adjust(x,method="BY")<p_thr))) / length(result_list)
sum(sapply(result_list,function(x)any(p.adjust(x,method="fdr")<p_thr))) / length(result_list)
