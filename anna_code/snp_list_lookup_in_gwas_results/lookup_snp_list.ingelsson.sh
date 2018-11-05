#!/bin/bash
python lookup_snp_list.py --var_list mhc/mhc.$2 --gwas_hits /oak/stanford/groups/euan/projects/ukbb/gwas/physical_activity_plink/v2/emmi/$1/$1.$2.continuous.assoc.linear --outf $1.chrom$2 --phenotype $1 --index_col 1
