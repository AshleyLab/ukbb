#!/bin/bash
python merge_and_annotate.py --exclude_list /oak/stanford/groups/euan/projects/ukbb/data/genetic_data/hrc/bad.snps.txt \
       --maf_prefix /oak/stanford/groups/euan/projects/ukbb/data/genetic_data/v2/EGAD00010001225/001/ukb_mfi_chr \
       --gwas_catalog /scratch/PI/euan/common/gwascatalog/gwas_catalog_may2017.slim.txt \
       --outf /oak/stanford/groups/euan/projects/ukbb/code/anna_code/postprocessing/top_snps_activity_merged.txt \
       --dir_to_merge /oak/stanford/groups/euan/projects/ukbb/code/anna_code/postprocessing/top_snps

