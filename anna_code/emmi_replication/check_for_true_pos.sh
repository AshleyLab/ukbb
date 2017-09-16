#!/bin/bash
#python check_for_true_pos.py --gs_file gs/chronotype_GS.csv \
#       --gs_chrom 1 \
#       --gs_pos 2 \
#       --ashley_lab_gwas_prefix /oak/stanford/groups/euan/projects/ukbb/gwas/physical_activity_plink/v2/emmi/chronotype/chronotype \
#       --ashley_lab_chrom 0 \
#       --ashley_lab_pos 2 \
#       --outf chronotype_tp_check.tsv


python check_for_true_pos.py --gs_file gs/grip_strength_GS.csv \
       --gs_chrom 1 \
       --gs_pos 2 \
       --ashley_lab_gwas_prefix /oak/stanford/groups/euan/projects/ukbb/gwas/physical_activity_plink/v2/emmi/grip_strength_to_bmi/grip_strength_to_bmi \
       --ashley_lab_chrom 0 \
       --ashley_lab_pos 2 \
       --outf grip_strength_to_bmi_tp_check.tsv



python check_for_true_pos.py --gs_file gs/grip_strength_GS.csv \
       --gs_chrom 1 \
       --gs_pos 2 \
       --ashley_lab_gwas_prefix /oak/stanford/groups/euan/projects/ukbb/gwas/physical_activity_plink/v2/emmi/grip_strength_to_weight/grip_strength_to_weight \
       --ashley_lab_chrom 0 \
       --ashley_lab_pos 2 \
       --outf grip_strength_to_weight_tp_check.tsv


python check_for_true_pos.py --gs_file gs/grip_strength_GS.csv \
       --gs_chrom 1 \
       --gs_pos 2 \
       --ashley_lab_gwas_prefix /oak/stanford/groups/euan/projects/ukbb/gwas/physical_activity_plink/v2/emmi/max_grip_strength/max_grip_strength \
       --ashley_lab_chrom 0 \
       --ashley_lab_pos 2 \
       --outf max_grip_strength_tp_check.tsv


python check_for_true_pos.py --gs_file gs/sedentary_GS.csv \
       --gs_chrom 1 \
       --gs_pos 2 \
       --ashley_lab_gwas_prefix /oak/stanford/groups/euan/projects/ukbb/gwas/physical_activity_plink/v2/emmi/rn_sedentary_behaviour/rn_sedentary_behaviour \
       --ashley_lab_chrom 0 \
       --ashley_lab_pos 2 \
       --outf rn_sedentary_behaviour_tp_check.tsv


python check_for_true_pos.py --gs_file gs/sleep_duration_GS.csv \
       --gs_chrom 1 \
       --gs_pos 2 \
       --ashley_lab_gwas_prefix /oak/stanford/groups/euan/projects/ukbb/gwas/physical_activity_plink/v2/emmi/sleep_duration/sleep_duration \
       --ashley_lab_chrom 0 \
       --ashley_lab_pos 2 \
       --outf sleep_duration_behaviour_tp_check.tsv


python check_for_true_pos.py --gs_file gs/METh_per_week_GS.csv \
       --gs_chrom 0 \
       --gs_pos 1 \
       --ashley_lab_gwas_prefix /oak/stanford/groups/euan/projects/ukbb/gwas/physical_activity_plink/v2/emmi/rn_total_METh_per_week/rn_total_METh_per_week \
       --ashley_lab_chrom 0 \
       --ashley_lab_pos 2 \
       --outf rn_total_METh_pwer_week_tp_check.tsv



