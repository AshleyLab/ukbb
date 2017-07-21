#!/bin/sh
#original covariates (sex & year of birth pulled from the ukb7454 table) 
#python aggregate_fields.py --tables /oak/stanford/groups/euan/projects/ukbb/data/bulk_data/ukb7454.tab \
#       --fields /oak/stanford/groups/euan/projects/ukbb/code/anna_code/gwas/ukb7454.fields \
#       --outf /oak/stanford/groups/euan/projects/ukbb/gwas/accelerometry_plink/v1/accelerometery_covariates

#adding in sex, year of birth, and principal components from population stratification -- EUROPEAN SUBJECTS 
python aggregate_fields.py --tables /oak/stanford/groups/euan/projects/ukbb/data/bulk_data/ukb7454.tab /oak/stanford/groups/euan/projects/ukbb/gwas/pop_strat/v2/euro/pca_results_v2_chrom1_euro.eigenvec \
       --fields /oak/stanford/groups/euan/projects/ukbb/code/anna_code/gwas/ukb7454.fields /oak/stanford/groups/euan/projects/ukbb/code/anna_code/gwas/pca.filtered.fields \
       --outf /oak/stanford/groups/euan/projects/ukbb/gwas/accelerometry_plink/v2/accelerometery_covariates
