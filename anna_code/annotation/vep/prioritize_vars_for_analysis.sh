#python prioritize_vars_for_analysis.py --vep_results activity_results/variant_effect_output.txt \
#       --top_snps /oak/stanford/groups/euan/projects/ukbb/code/anna_code/postprocessing/self_reported_vs_measured/top_snps_activity_merged.txt \
#       --impact_column 4 \
#       --outf prioritized.activity.vars

python prioritize_vars_for_analysis.py --vep_results fitness_results/variant_effect_output.txt \
       --top_snps /oak/stanford/groups/euan/projects/ukbb/code/anna_code/postprocessing/fitness_hits.tsv \
       --impact_column 4 \
       --outf prioritized.fitness.vars
