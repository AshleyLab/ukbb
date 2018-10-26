#!/bin/bash
python lookup_snp_list.py --var_list mhc/mhc.$2 --gwas_hits ../../../gwas/physical_activity_plink/v2/results_imputed/X$1/$1.$2.continuous.$1.glm.linear --outf $1.chrom$2 --phenotype $1
