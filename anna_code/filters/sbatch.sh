for chrom in `seq 1 22`
do

    sbatch -J "called_vs_imp.$chrom" -o logs/$chrom.o -e logs/$chrom.e -p euan,owners panel_intersect.sh $chrom
done
