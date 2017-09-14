#!/bin/bash
#2/3 discovery, 1/3 replication for maf thresholds:

# maf < 0.005 
#python validate_discovery_and_replication.py --discovery_dir /oak/stanford/groups/euan/projects/ukbb/gwas/physical_activity_plink/v2/discovery \
#       --replication_dir /oak/stanford/groups/euan/projects/ukbb/gwas/physical_activity_plink/v2/validation \
#       --outf discovery.replication.maf_lt_0.005 \
#       --bidirectional \
#       --features DurationWalkingForPleasure OverallAccelerationAverage TimeSpentOutdoorsSummer  \
#       --maf_thresh \
#       --maf_min 0.001 \
#       --maf_max 0.005 \
#       --maf_metadata /oak/stanford/groups/euan/projects/ukbb/data/genetic_data/v2/EGAD00010001225/001 
       

# maf in range [0.005,0.01] 
python validate_discovery_and_replication.py --discovery_dir /oak/stanford/groups/euan/projects/ukbb/gwas/physical_activity_plink/v2/discovery \
       --replication_dir /oak/stanford/groups/euan/projects/ukbb/gwas/physical_activity_plink/v2/validation \
       --outf discovery.replication.maf.range.0.005.to.0.01 \
       --bidirectional \
       --features DurationWalkingForPleasure OverallAccelerationAverage TimeSpentOutdoorsSummer \
       --maf_thresh \
       --maf_min 0.005 \
       --maf_max 0.01 \
       --maf_metadata /oak/stanford/groups/euan/projects/ukbb/data/genetic_data/v2/EGAD00010001225/001 


       
# maf > 0.01 
#python validate_discovery_and_replication.py --discovery_dir /oak/stanford/groups/euan/projects/ukbb/gwas/physical_activity_plink/v2/discovery \
#       --replication_dir /oak/stanford/groups/euan/projects/ukbb/gwas/physical_activity_plink/v2/validation \
#       --outf discovery.replication.maf_gt_0.01 \
#       --bidirectional \
#       --features DurationWalkingForPleasure OverallAccelerationAverage TimeSpentOutdoorsSummer \
#       --maf_min 0.01 \
#       --maf_max 0.5 \
#       --maf_metadata /oak/stanford/groups/euan/projects/ukbb/data/genetic_data/v2/EGAD00010001225/001
