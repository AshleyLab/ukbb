#!/bin/sh
python extract_cluster_trajectories.py --base_dir /scratch/PI/euan/projects/ukbb/data/accelerometry/aligned\
       --subject_to_cluster /scratch/PI/euan/projects/ukbb/data/accelerometry/aligned/cluster_mapping.tsv\
       --cluster_of_interest $1\
       --outf /scratch/PI/euan/projects/ukbb/data/accelerometry/aligned/mean_trajectories/$1.tsv
