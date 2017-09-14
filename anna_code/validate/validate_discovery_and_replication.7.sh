#!/bin/bash
#compare at locus level, w/ locus being 500 kb flank (250 kb upstream, 250 kb downstream)
python validate_discovery_and_replication.py --discovery_dir /oak/stanford/groups/euan/projects/ukbb/gwas/physical_activity_plink/v2/discovery \
       --replication_dir /oak/stanford/groups/euan/projects/ukbb/gwas/physical_activity_plink/v2/validation \
       --outf locus.250 \
       --bidirectional \
       --features DurationWalkingForPleasure OverallAccelerationAverage TimeSpentOutdoorsSummer \
       --make_flank_bed \
       --flank_size 250000
       
#compare at locus level, w/ locus being 1000 kb flank (500 kb upstream, 500 kb downstream)
python validate_discovery_and_replication.py --discovery_dir /oak/stanford/groups/euan/projects/ukbb/gwas/physical_activity_plink/v2/discovery \
       --replication_dir /oak/stanford/groups/euan/projects/ukbb/gwas/physical_activity_plink/v2/validation \
       --outf locus.500 \
       --bidirectional \
       --features DurationWalkingForPleasure OverallAccelerationAverage TimeSpentOutdoorsSummer \
       --make_flank_bed \
       --flank_size 500000
       
