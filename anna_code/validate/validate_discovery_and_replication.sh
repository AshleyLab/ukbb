#!/bin/bash

#2/3 discovery, 1/3 replication 
#python validate_discovery_and_replication.py --discovery_dir /oak/stanford/groups/euan/projects/ukbb/gwas/physical_activity_plink/v2/discovery \
#       --replication_dir /oak/stanford/groups/euan/projects/ukbb/gwas/physical_activity_plink/v2/validation \
#       --outf discovery.replication \
#       --bidirectional

#extra covariates of BMI, height, and batch 
#python validate_discovery_and_replication.py --discovery_dir /oak/stanford/groups/euan/projects/ukbb/gwas/physical_activity_plink/v2/results_extracovars \
#       --replication_dir /oak/stanford/groups/euan/projects/ukbb/gwas/physical_activity_plink/v2/results \
#       --outf extracovars \
#       --bidirectional 


#extra covariates of BMI, height, Batch, 10 pc's 
#python validate_discovery_and_replication.py --discovery_dir /oak/stanford/groups/euan/projects/ukbb/gwas/physical_activity_plink/v2/results_extracovars_10pcs \
#       --replication_dir /oak/stanford/groups/euan/projects/ukbb/gwas/physical_activity_plink/v2/results \
#       --outf 10pcs \
#       --bidirectional 

#not using a discovery set for 23&Me features
python validate_discovery_and_replication.py --discovery_dir /oak/stanford/groups/euan/projects/ukbb/gwas/physical_activity_plink/v2/results \
       --outf for23andme 
