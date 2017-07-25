#!/bin/bash
geno=/oak/stanford/groups/euan/projects/ukbb/data/genetic_data/v2/plink/ukb_imp_chr10_v2
pheno="/oak/stanford/groups/euan/projects/ukbb/da_dh/pheno_matrices/Predicted_HR_at_WD=100_simple_euro.txt"
phenoname=Residuals_HR_pred_100
bolt --bfile=$geno --phenoFile=$pheno --phenoCol=$phenoname --reml --remlNoRefine --numThreads=4 > ${phenoname}.out
