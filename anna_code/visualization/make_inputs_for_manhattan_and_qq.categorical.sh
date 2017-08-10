#!/bin/bash

python make_inputs_for_manhattan_and_qq.py --exclude_pickle /oak/stanford/groups/euan/projects/ukbb/code/anna_code/useful_pickles/exclude_dict.p \
       --result_dir /oak/stanford/groups/euan/projects/ukbb/gwas/accelerometry_plink/v2/results/$1 \
       --tokeep SNP CHR BP P \
       --ending .categorical.qassoc \
       --outf $1.aggregate.maf.0.01 \
       --maf_filter 0.01 \
       --maf_pickle_prefix /oak/stanford/groups/euan/projects/ukbb/code/anna_code/useful_pickles/maf_

python make_inputs_for_manhattan_and_qq.py --exclude_pickle /oak/stanford/groups/euan/projects/ukbb/code/anna_code/useful_pickles/exclude_dict.p \
       --result_dir /oak/stanford/groups/euan/projects/ukbb/gwas/accelerometry_plink/v2/results/$1 \
       --tokeep SNP CHR BP P \
       --ending .categorical.qassoc \
       --outf $1.aggregate \
