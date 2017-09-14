#!/bin/bash

#2/3 discovery, 1/3 replication 
python validate_discovery_and_replication.py --discovery_dir /oak/stanford/groups/euan/projects/ukbb/gwas/physical_activity_plink/v2/discovery \
       --replication_dir /oak/stanford/groups/euan/projects/ukbb/gwas/physical_activity_plink/v2/validation \
       --outf discovery.replication \
       --bidirectional
