#!/bin/bash
#2/3 discovery, 1/3 replication for directly genotyped SNPs only 
python validate_discovery_and_replication.py --discovery_dir /oak/stanford/groups/euan/projects/ukbb/gwas/physical_activity_plink/v2/discovery \
       --replication_dir /oak/stanford/groups/euan/projects/ukbb/gwas/physical_activity_plink/v2/validation \
       --outf discovery.replication.calledonly \
       --bidirectional \
       --directly_genotyped \
       --directly_genotyped_metadata /oak/stanford/groups/euan/projects/ukbb/data/genetic_data/v2/EGAD00010001226/001 \
       --features DurationWalkingForPleasure OverallAccelerationAverage TimeSpentOutdoorsSummer 
       
