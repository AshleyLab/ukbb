#!/bin/sh
python aggregate_fields.py --tables /scratch/PI/euan/projects/ukbb/data/bulk_data/ukb7454.tab \
       --fields /scratch/PI/euan/projects/ukbb/code/anna_code/gwas/ukb7454.fields \
       --outf /scratch/PI/euan/projects/ukbb/gwas/accelerometry_plink/v1/accelerometery_covariates
