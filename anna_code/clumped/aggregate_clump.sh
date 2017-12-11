#!/bin/bash
python aggregate_clump.py --sourcef /oak/stanford/groups/euan/projects/ukbb/gwas/physical_activity_plink/v2/time_outdoors/$1/all.clumped \
       --outf $HOME/$1.clumped \
       --verify_parts /oak/stanford/groups/euan/projects/ukbb/gwas/physical_activity_plink/v2/time_outdoors/$1/$1. .continuous.$1.glm.linear

