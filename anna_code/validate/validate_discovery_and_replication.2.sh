#!/bin/bash

#extra covariates of BMI, height, and batch 
python validate_discovery_and_replication.py --discovery_dir /oak/stanford/groups/euan/projects/ukbb/gwas/physical_activity_plink/v2/results \
       --replication_dir /oak/stanford/groups/euan/projects/ukbb/gwas/physical_activity_plink/v2/results_extracovars \
       --features DurationWalkingForPleasure OverallAccelerationAverage TimeSpentOutdoorsSummer \
       --pval_validation_thresh 5e-08 \
       --outf extracovars \
       --bidirectional 
