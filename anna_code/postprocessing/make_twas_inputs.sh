#!/bin/bash 
python make_twas_inputs.py --frequencies /oak/stanford/groups/euan/projects/ukbb/data/genetic_data/v2/EGAD00010001225/001/ukb_mfi_chr$1\_v2.txt \
       --gwas_output da_dh_jul21_gwas.chr$1.txt --keep_columns SNP A1 A2 STAT --exclude_list /oak/stanford/groups/euan/projects/ukbb/data/genetic_data/hrc/bad.snps.txt --outdir /oak/stanford/groups/euan/projects/ukbb/twas/twas_inputs_da_dh_jul21


