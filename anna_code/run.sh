#!/bin/bash
#gconv $PI_SCRATCH/projects/ukbb/chrom10.cal c10 ped -r
#plink --file $PI_SCRATCH/projects/ukbb/data/genetic_data/calls/plink/chrom10 --make-bed --out $PI_SCRATCH/projects/ukbb/data/genetic_data/calls/plink/chrom10
plink --file $PI_SCRATCH/projects/ukbb/data/genetic_data/calls/plink/chrom10 --make-bed --out chrom10

