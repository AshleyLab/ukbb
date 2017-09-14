#!/bin/bash

#not using a discovery set for 23&Me features
python validate_discovery_and_replication.py --discovery_dir /oak/stanford/groups/euan/projects/ukbb/gwas/physical_activity_plink/v2/results \
       --outf for23andme 
