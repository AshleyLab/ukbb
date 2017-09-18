for i in `cat fields.categorical`
do
    for chrom in `seq 1 22`
    do
	sbatch -J "$i.$chrom.categorical" -o logs2/$i.$chrom.o -e logs2/$i.$chrom.e -p akundaje,owners  --time=08:00:00 --mem=20000  plink_script.imputed.categorical.sh $i $chrom
    done
done
