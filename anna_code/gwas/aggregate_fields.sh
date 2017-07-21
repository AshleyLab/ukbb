#!/bin/sh

#ACTIVITY FIELDS 
#python aggregate_fields.py --tables /scratch/PI/euan/projects/ukbb/data/bulk_data/ukb8524.tab /scratch/PI/euan/projects/ukbb/data/bulk_data/ukb7454.tab /scratch/PI/euan/projects/ukbb/gwas/accelerometry_plink/v1/transitions/transitions.10.tsv /scratch/PI/euan/projects/ukbb/gwas/accelerometry_plink/v1/transitions/transitions.25.tsv /scratch/PI/euan/projects/ukbb/gwas/accelerometry_plink/v1/dwt/dwt_features_small.txt \
#       --fields /scratch/PI/euan/projects/ukbb/code/anna_code/gwas/ukb8524.fields /scratch/PI/euan/projects/ukbb/code/anna_code/gwas/ukb7454.fields /scratch/PI/euan/projects/ukbb/code/anna_code/gwas/transitions.10.fields /scratch/PI/euan/projects/ukbb/code/anna_code/gwas/transitions.25.fields /scratch/PI/euan/projects/ukbb/code/anna_code/gwas/dwt_features_small.fields \
#       --outf /scratch/PI/euan/projects/ukbb/gwas/accelerometry_plink/v1/accelerometery_aggregate_phenotypes


#ETHNIC BACKGROUND

python aggregate_fields.py --tables  /scratch/PI/euan/projects/ukbb/data/bulk_data/ukb7454.tab  \
       --fields /scratch/PI/euan/projects/ukbb/code/anna_code/gwas/ethnicity.fields  \
       --outf /scratch/PI/euan/projects/ukbb/gwas/accelerometry_plink/v1/ethnicity_aggregate_phenotypes

