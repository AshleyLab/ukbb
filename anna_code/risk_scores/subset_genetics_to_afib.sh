#!/bin/bash
#plink2 --bfile /oak/stanford/groups/euan/projects/ukbb/data/genetic_data/v2/plink_small/ukb_imp_chr$1\_v2 --keep afib.cohort.training.txt --make-bed --out afib.cohort.train.$1
plink2 --bfile /oak/stanford/groups/euan/projects/ukbb/data/genetic_data/v2/plink_small/ukb_imp_chr$1\_v2 --keep afib.cohort.test.txt  --make-bed --out afib.cohort.test.$1


