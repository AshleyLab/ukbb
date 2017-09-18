#for i in `cat fields.continuous`
for i in 0_1 10_11 11_12 1_2 12_13 13_14 14_15 15_16 16_17 17_18 18_19 19_20 2000mg 20_21 21_22 22_23 2_3 23_24 3_4 4_5 5_6 6_7 7_8 8_9 9_10 9mg 
do
    for chrom in `seq 1 22`
    do
	sbatch -J "$i.$chrom.continuous" -o logs2/$i.$chrom.o -e logs2/$i.$chrom.e -p akundaje,owners  --time=08:00:00 --mem=20000  plink_script.imputed.continuous.residuals_qnorm.sh $i $chrom
    done
done
