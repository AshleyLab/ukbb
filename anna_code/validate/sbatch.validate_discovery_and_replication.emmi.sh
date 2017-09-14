for feature in chronotype grip_strength_to_bmi grip_strength_to_weight max_grip_strength rn_sedentary_behaviour rn_total_METh_per_week sleep_duration
do    
    sbatch -J "filter.$feature" -o logs/$feature.o -e logs/$feature.e -p euan,owners  validate_discovery_and_replication.emmi.sh $feature 
done
