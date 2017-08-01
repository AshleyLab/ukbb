#!/bin/bash 
python make_inputs_for_manhattan_and_qq.py --exclude_list /oak/stanford/groups/euan/projects/ukbb/data/genetic_data/hrc/bad.snps.txt \
       --result_dir /oak/stanford/groups/euan/projects/ukbb/gwas/accelerometry_plink/v2/results/$1 \
       --tokeep SNP CHR BP P \
       --ending .categorical.qassoc \
       --outf $1.aggregate
