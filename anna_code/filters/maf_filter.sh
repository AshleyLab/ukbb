#filter the rn_METh PLINK hits for Emmi's phenotype verification
python maf_filter.py --input /oak/stanford/groups/euan/projects/ukbb/gwas/physical_activity_plink/v2/emmi/rn_total_METh_per_week/all.assoc.linear \
       --thresh 5e-08 \
       --outf /oak/stanford/groups/euan/projects/ukbb/gwas/physical_activity_plink/v2/emmi/rn_total_METh_per_week/all.assoc.linear.filtered

