#!/bin/bash
#extra covariates of BMI, height, Batch, 10 pc's 
python validate_discovery_and_replication.py --discovery_dir /oak/stanford/groups/euan/projects/ukbb/gwas/physical_activity_plink/v2/results \
       --replication_dir /oak/stanford/groups/euan/projects/ukbb/gwas/physical_activity_plink/v2/results_extracovars_10pcs \
       --pval_validation_thresh 5e-08 \
       --features DurationWalkingForPleasure OverallAccelerationAverage TimeSpentOutdoorsSummer \
       --outf 10pcs \
       --bidirectional 

