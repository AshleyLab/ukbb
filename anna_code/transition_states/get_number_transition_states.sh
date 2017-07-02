#!/bin/sh
start_index=$1
end_index=$2
python get_number_transition_states.py --subject_file /scratch/PI/euan/projects/ukbb/data/accelerometry/aligned/aligned.txt\
       --mg_cutoff 25\
       --outf /scratch/PI/euan/projects/ukbb/gwas/accelerometry_plink/v1/transition_states.25.$start_index.tsv\
       --start_index $start_index\
       --end_index $end_index

