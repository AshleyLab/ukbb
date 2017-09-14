python add_effect_size_info.py --snp_hit_file /oak/stanford/groups/euan/projects/ukbb/code/anna_code/postprocessing/top_snps_activity_merged.txt \
       --effect_col_name MeanEffectSize \
       --feature_col_name Features \
       --file_to_annotate prioritized.activity.vars \
       --outf prioritized.activity.vars.with.effect

python add_effect_size_info.py --snp_hit_file /oak/stanford/groups/euan/projects/ukbb/code/anna_code/postprocessing/fitness_hits.tsv \
       --effect_col_name effect \
       --feature_col_name trait \
       --file_to_annotate prioritized.fitness.vars \
       --outf prioritized.fitness.vars.with.effect

