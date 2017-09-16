echo "chronotype"
python collapse_clumps.py --hp_dir /oak/stanford/groups/euan/projects/ukbb/gwas/physical_activity_plink/v2/emmi/chronotype/chronotype \
       --gs_snps gs/chronotype.snp \
       --called_snps plink_ld_filter/emmi.chronotype.filtered
echo "sleep duration"
python collapse_clumps.py --hp_dir /oak/stanford/groups/euan/projects/ukbb/gwas/physical_activity_plink/v2/emmi/sleep_duration/sleep_duration \
       --gs_snps gs/sleep.snp \
       --called_snps plink_ld_filter/emmi.sleep_duration.filtered

echo "grip strength to weight"
python collapse_clumps.py --hp_dir /oak/stanford/groups/euan/projects/ukbb/gwas/physical_activity_plink/v2/emmi/grip_strength_to_weight \
       --gs_snps gs/grip.snp \
       --called_snps plink_ld_filter/emmi.grip_strength_to_weight.filtered

echo "sedentary behavior"
python collapse_clumps.py --hp_dir /oak/stanford/groups/euan/projects/ukbb/gwas/physical_activity_plink/v2/emmi/rn_sedentary_behaviour/rn_sedentary_behaviour \
       --gs_snps gs/sed.snp \
       --called_snps plink_ld_filter/emmi.rn_sedentary_behaviour.filtered




#chronotype.snp  grip.snp  sed.snp  sleep.snp

