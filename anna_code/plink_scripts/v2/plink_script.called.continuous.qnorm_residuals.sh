#!/bin/sh
mkdir /oak/stanford/groups/euan/projects/ukbb/gwas/physical_activity_plink/qc/called_qnorm_residuals/$1 
$HOME/apps/plink/plink2 --bfile /oak/stanford/groups/euan/projects/ukbb/data/genetic_data/v2/plink_dir_genotype/ukb_cal_chr$2\_v2 --pheno overall_acceleration_average_qnorm_residuals.txt   --input-missing-phenotype -1000  --out /oak/stanford/groups/euan/projects/ukbb/gwas/physical_activity_plink/qc/called_qnorm_residuals/$1/$1.$2.continuous  --pheno-name $1 --keep euro_minus_exclusion_minus_firstdegree.txt  --glm
