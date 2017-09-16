base_dir=/oak/stanford/groups/euan/projects/ukbb/code/anna_code/emmi_replication/pval_filter
#for feature in emmi.chronotype  emmi.grip_strength_to_bmi  emmi.grip_strength_to_weight  emmi.max_grip_strength  emmi.rn_sedentary_behaviour  emmi.rn_total_METh_per_week  emmi.sleep_duration
for feature in emmi.sleep_duration
do
    full_path=$base_dir/$feature
    sbatch -J "badsnps.$feature" -o logs/$feature.o -e logs/$feature.e -p euan,owners exclude_bad_snps.sh $full_path
done

	       
