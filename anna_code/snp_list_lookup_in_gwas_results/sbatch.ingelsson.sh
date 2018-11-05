for phenotype in chronotype grip_strength_to_weight rn_sedentary_behaviour sleep_duration grip_strength_to_bmi max_grip_strength rn_total_METh_per_week
do
    for chrom in `seq 1 22`
    do
        sbatch -J "$phenotype.$chrom.continuous" -o logs/$phenotype.$chrom.o -e logs/$phenotype.$chrom.e -p akundaje,euan,owners  lookup_snp_list.ingelsson.sh $phenotype $chrom
	#sbatch -J "$phenotype.$chrom.continuous" -o logs/$phenotype.$chrom.o -e logs/$phenotype.$chrom.e -p akundaje,euan,owners  lookup_snp_list.categorical.sh $phenotype $chrom
    done
done
