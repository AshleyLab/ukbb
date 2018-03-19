python mr_input_df_from_plink.py --exposure_files /oak/stanford/groups/euan/projects/ukbb/code/anna_code/mr/mr_snp_subsets/* \
    --outcomes DeathStatus DeathAge CvdStatus AssessedAge AliveAt70 AliveAt65 \
    --exposure_prefix /oak/stanford/groups/euan/projects/ukbb/ukbb/gwas/physical_activity_plink/v2/results_imputed \
    --outcome_prefix /oak/stanford/groups/euan/projects/ukbb/ukbb/gwas/physical_activity_plink/v2/mr_outcomes \
    --outf gwas.summary.test
