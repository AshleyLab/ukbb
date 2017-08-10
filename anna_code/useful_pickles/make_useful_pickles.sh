#!/bin/bash
python make_useful_pickles.py --maf_base /oak/stanford/groups/euan/projects/ukbb/data/genetic_data/v2/EGAD00010001225/001/ukb_mfi_chr  --chrom $1 --outf maf_$1.p

