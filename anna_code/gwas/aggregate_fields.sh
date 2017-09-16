#!/bin/sh

#QC -- OVERALL ACCELERATION AVERAGE 

#ACTIVITY FIELDS 
python aggregate_fields.py --tables /oak/stanford/groups/euan/projects/ukbb/data/bulk_data/ukb8524.tab /oak/stanford/groups/euan/projects/ukbb/data/bulk_data/ukb7454.tab /oak/stanford/groups/euan/projects/ukbb/code/anna_code/transition_states/transitions/transitions.10.tsv /oak/stanford/groups/euan/projects/ukbb/code/anna_code/transition_states/transitions/transitions.25.tsv /oak/stanford/groups/euan/projects/ukbb/code/anna_code/dwt/dwt_features_small.txt \
       --fields /oak/stanford/groups/euan/projects/ukbb/code/anna_code/gwas/ukb8524.fields /oak/stanford/groups/euan/projects/ukbb/code/anna_code/gwas/ukb7454.fields /oak/stanford/groups/euan/projects/ukbb/code/anna_code/gwas/transitions.10.fields /oak/stanford/groups/euan/projects/ukbb/code/anna_code/gwas/transitions.25.fields /oak/stanford/groups/euan/projects/ukbb/code/anna_code/gwas/dwt_features_small.fields \
       --outf accelerometery_aggregate_phenotypes


#ETHNIC BACKGROUND
#python aggregate_fields.py --tables  /oak/stanford/groups/euan/projects/ukbb/data/bulk_data/ukb7454.tab  \
#       --fields /oak/stanford/groups/euan/projects/ukbb/code/anna_code/gwas/ethnicity.fields  \
#       --outf /oak/stanford/groups/euan/projects/ukbb/gwas/accelerometry_plink/v1/ethnicity_aggregate_phenotypes


#DURATION OF VIGOROUS ACTIVITY
#python aggregate_fields.py --tables /oak/stanford/groups/euan/projects/ukbb/data/bulk_data/ukb7454.tab \
#       --fields /oak/stanford/groups/euan/projects/ukbb/code/anna_code/gwas/duration.vigorous.activity.fields \
#       --outf accelerometery_aggregate_phenotypes.DurationVigorous.Activity.txt



