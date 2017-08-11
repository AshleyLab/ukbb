python merge_motifs_with_vep.py --motif_hits activity__1_1e-06/motif_hits_collapsed.txt \
       --vep_annotation /oak/stanford/groups/euan/projects/ukbb/code/anna_code/annotation/vep/activity_results/variant_effect_output.txt \
       --outf activity_regulatory_annotation.txt


python merge_motifs_with_vep.py --motif_hits fitness__1_1e-06/motif_hits_collapsed.txt \
       --vep_annotation /oak/stanford/groups/euan/projects/ukbb/code/anna_code/annotation/vep/fitness_results/variant_effect_output.txt \
       --outf fitness_regulatory_annotation.txt

