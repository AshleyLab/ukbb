#for chrom in `seq 1 22` 
#do
#    sbatch -J "$chrom.test" -o logs/$chrom.test.o -e logs/$chrom.test.e -p euan,owners  --time=48:00:00 --mem=20000  subset_genetics_to_afib.sh $chrom
#done

