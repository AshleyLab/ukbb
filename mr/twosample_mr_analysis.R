library('MendelianRandomization')

get_twosample_mr_input<-function(effects_data){
  x_inds = effects_data$Source!="RHR"
  bx = -effects_data[x_inds,"Effect"]
  bxse = effects_data[x_inds,"SE"]
  y_inds = effects_data$Source=="RHR"
  by = -effects_data[y_inds,"Effect"]
  byse = effects_data[y_inds,"SE"]
  all(effects_data$SNP[x_inds] == effects_data$SNP[y_inds])
  snp_names = effects_data$SNP[x_inds]
  names(bx) = snp_names;names(bxse) = snp_names
  names(by) = snp_names; names(byse) = snp_names
  return(data.frame(bx=bx,bxse=bxse,by=by,byse=byse))
}

# recovery analysis
input_file = '/Users/David/Desktop/ukbb/mr/gwas_effects/recovery_snps_all_effects.txt'
effects_data = read.delim(input_file)
mr_data = get_twosample_mr_input(effects_data)
mr_in = mr_input(mr_data$bx,mr_data$bxse,mr_data$by,mr_data$byse)
mr_plot(mr_in)
mr_allmethods(mr_in)
mr_in = mr_input(abs(mr_data$bx),mr_data$bxse,abs(mr_data$by),mr_data$byse)
mr_plot(mr_in)
mr_allmethods(mr_in)

input_file = '/Users/David/Desktop/ukbb/mr/gwas_effects/exercise_snps_all_effects.txt'
effects_data = read.delim(input_file)
mr_data = get_twosample_mr_input(effects_data)
mr_in = mr_input(mr_data$bx,mr_data$bxse,mr_data$by,mr_data$byse)
mr_plot(mr_in)
mr_allmethods(mr_in)

# analyze the exercise snps bu their effects
low_var_snps = rownames(mr_data)[(mr_data$bxse<0.01)]
high_var_snps = rownames(mr_data)[(mr_data$bxse>=0.01)]

# interpret the results using FUMA's output
fuma_path = "Desktop/ukbb/gwas_interpretation/fuma/pooled_p1_5e8_p2_1e4_ld_0.6_with_maf/"
get_fuma_gwascatalog_results<-function(path){
  res = read.delim(paste(path,"gwascatalog.txt",sep=''),stringsAsFactors = F)
  res = res[,c("GenomicLocus","IndSigSNP","snp","Trait","Region","ReportedGene","MappedGene","SNPs","P")]
  return(res)
}
get_fuma_gene_sets_results<-function(path){
  res = read.delim(paste(path,"GS.txt",sep=''),stringsAsFactors = F)
  return(res)
}

gwastcat_res = get_fuma_gwascatalog_results(fuma_path)
trait2leadsnps = tapply(gwastcat_res$IndSigSNP,gwastcat_res$Trait,unique)
trait2leadsnps[["NA"]] = setdiff(rownames(mr_data),gwastcat_res$IndSigSNP)
high_var_snps_gwascat = gwastcat_res[is.element(gwastcat_res$IndSigSNP,set = high_var_snps),]
low_var_snps_gwascat = gwastcat_res[is.element(gwastcat_res$IndSigSNP,set = low_var_snps),]
sort(table(low_var_snps_gwascat$Trait))

# MR using the resting HR snps
curr_snps = trait2leadsnps$`Resting heart rate`
inds = is.element(rownames(mr_data),set=curr_snps)

par(mfrow=c(1,2))
boxplot(mr_data$bxse~inds, ylab = "Effect SE, Exercise HR", 
        col = c("darkred","cyan"), names = c("Other","RHR SNPs"))
boxplot(abs(mr_data$by)~inds, ylab = "Effect size, RHR",
        col = c("darkred","cyan"), names = c("Other","RHR SNPs"))

fisher.test(table(inds,mr_data$bxse<0.01))

mr_in = mr_input(mr_data$bx[inds],mr_data$bxse[inds],mr_data$by[inds],mr_data$byse[inds])
mr_plot(mr_in)
mr_allmethods(mr_in)

# SANITY CHECK
# Comparison to the old, direct analysis
# For this analysis to work, first run the non-activity
# analysis in the direct_mr_analysis.R script
# What is needed from that script is the genotype and the linear
# regression analyses.
direct_analysis_in = mr_analysis_res[["fitness_hr,fitness_hr,rhr"]]$multivar_input
shared_snps = intersect(names(slot(direct_analysis_in,"betaX")),snp_names)
bx1 = bx[shared_snps]
bx2 = slot(direct_analysis_in,"betaX")[shared_snps]
bxse1 = bxse[shared_snps]
bxse2 = slot(direct_analysis_in,"betaXse")[shared_snps]



