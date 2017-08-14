#!/bin/bash
python panel_intersect.py --imputed_snps /oak/stanford/groups/euan/projects/ukbb/data/genetic_data/v2/EGAD00010001225/001/ukb_mfi_chr$1\_v2.txt \
       --called_snps /oak/stanford/groups/euan/projects/ukbb/data/genetic_data/v2/EGAD00010001226/001/ukb_snp_chr$1\_v2.bim \
       --outf $1 
