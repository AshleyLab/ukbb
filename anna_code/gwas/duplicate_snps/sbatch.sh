for chrom in `seq 1 22`
do
    sbatch -J "duplicate_snps$chrom" -o logs/duplicate.$chrom.o -e logs/duplicate.$chrom.e -p euan,owners get_list_of_duplicate_snps.sh $chrom
done

