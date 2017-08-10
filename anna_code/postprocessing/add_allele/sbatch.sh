for chrom in `seq 5 21`
do
    sbatch -J "add_allele_info.$chrom" -o logs/add_allele_info.$chrom.o -e logs/add_allele_info.$chrom.e -p euan,owners add_allele_to_gwas_results.sh $chrom
done
