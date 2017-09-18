for i in OverallAccelerationAverage
do
    for chrom in `seq 1 22`
    do
	sbatch -J "$i.$chrom.continuous" -o logs2/$i.$chrom.o -e logs2/$i.$chrom.e -p euan  --time=08:00:00  plink_script.called.continuous.residuals_qnorm.sh $i $chrom
    done
done
