for i in OverallAccelerationAverage
do
    for chrom in `seq 1 22`
    do
	sbatch -J "$i.$chrom.continuous" -o logs4/$i.$chrom.o -e logs4/$i.$chrom.e -p akundaje,owners  --time=08:00:00 --mem=20000  plink_script.imputed.continuous.sh $i $chrom
    done
done
