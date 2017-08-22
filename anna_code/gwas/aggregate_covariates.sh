#!/bin/sh
#original covariates (sex & year of birth pulled from the ukb7454 table) 
#python aggregate_fields.py --tables /oak/stanford/groups/euan/projects/ukbb/data/bulk_data/ukb7454.tab \
#       --fields /oak/stanford/groups/euan/projects/ukbb/code/anna_code/gwas/ukb7454.fields \
#       --outf /oak/stanford/groups/euan/projects/ukbb/gwas/accelerometry_plink/v1/accelerometery_covariates

#adding in sex, year of birth, and principal components from population stratification -- EUROPEAN SUBJECTS 
#python aggregate_fields.py --tables /oak/stanford/groups/euan/projects/ukbb/data/bulk_data/ukb7454.tab /oak/stanford/groups/euan/projects/ukbb/gwas/pop_strat/v2/euro/pca_results_v2_chrom1_euro.eigenvec batches.recoded \
#       --fields /oak/stanford/groups/euan/projects/ukbb/code/anna_code/gwas/ukb7454.fields /oak/stanford/groups/euan/projects/ukbb/code/anna_code/gwas/pca.filtered.fields batch.fields \
#       --outf covariates.txt


#adding in sex,  year of birth, first 5 pc's, standing height, and BMI
#python aggregate_fields.py --tables /oak/stanford/groups/euan/projects/ukbb/data/bulk_data/ukb7454.tab /oak/stanford/groups/euan/projects/ukbb/data/bulk_data/ukb7917.tab /oak/stanford/groups/euan/projects/ukbb/gwas/pop_strat/v2/euro/pca_results_v2_chrom1_euro.10pcs.eigenvec \
#       --fields /oak/stanford/groups/euan/projects/ukbb/code/anna_code/gwas/ukb7454.age.sex.fields /oak/stanford/groups/euan/projects/ukbb/code/anna_code/gwas/ukb7917.height.bmi.fields /oak/stanford/groups/euan/projects/ukbb/code/anna_code/gwas/pca.filtered.10pcs.fields \
#       --outf accelerometry_covariates.augmented.txt


#adding in sex,  year of birth, first 5 pc's, standing height, and BMI
python aggregate_fields.py --tables /oak/stanford/groups/euan/projects/ukbb/data/bulk_data/ukb7454.tab /oak/stanford/groups/euan/projects/ukbb/data/bulk_data/ukb7917.tab /oak/stanford/groups/euan/projects/ukbb/gwas/pop_strat/v2/euro/10pcs/pca_results_v2_chrom1_10pcseuro.eigenvec batches.recoded \
       --fields /oak/stanford/groups/euan/projects/ukbb/code/anna_code/gwas/ukb7454.age.sex.fields /oak/stanford/groups/euan/projects/ukbb/code/anna_code/gwas/ukb7917.height.bmi.fields /oak/stanford/groups/euan/projects/ukbb/code/anna_code/gwas/pca.filtered.fields batch.fields \
       --outf covariates.augmented.txt
