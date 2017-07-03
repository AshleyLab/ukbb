#!/bin/sh
module load R/3.2.0 
for entry in `cat /scratch/PI/euan/projects/ukbb/data/accelerometry/aligned/parts/x$1`
do
    Rscript extract_features_small.R $entry /scratch/PI/euan/projects/ukbb/gwas/accelerometry_plink/v1/dwt_features_small.$1.txt
    echo $entry 
done


