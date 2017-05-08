#!/bin/sh
module load R
Rscript feature_extract_accelerometry.R /scratch/PI/euan/projects/ukbb/data/accelerometry/aligned/x$1
