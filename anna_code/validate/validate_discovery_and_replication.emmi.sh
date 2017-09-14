#!/bin/bash

#Emmi's phenotypes -- no replication set for this analysis. 
python validate_discovery_and_replication.py --discovery_dir /oak/stanford/groups/euan/projects/ukbb/gwas/physical_activity_plink/v2/emmi \
       --outf emmi \
       --features $1

