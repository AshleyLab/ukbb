#!/bin/sh
mkdir /oak/stanford/groups/euan/projects/ukbb/gwas/physical_activity_plink/qc/imputed_residuals_qnorm/$1 
$HOME/apps/plink/plink2 --bfile /oak/stanford/groups/euan/projects/ukbb/data/genetic_data/v2/plink_small/ukb_imp_chr$2\_v2  --pheno overall_acceleration_average_residuals_qnorm.txt --input-missing-phenotype -1000 --out /oak/stanford/groups/euan/projects/ukbb/gwas/physical_activity_plink/qc/imputed_residuals_qnorm/$1/$1.$2.continuous  --pheno-name $1 --keep euro_minus_exclusion_minus_firstdegree.txt --glm
