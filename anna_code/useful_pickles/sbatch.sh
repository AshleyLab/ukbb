for chrom in `seq 4 5`
do
    sbatch -J "maf_pickle$chrom" -o logs/maf_pickle.$chrom.o -e logs/maf_pickle.$chrom.e -p euan --mem=10000 --time=04:00:00 make_useful_pickles.sh $chrom
done
