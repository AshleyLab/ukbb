for chrom in `seq 1 22` X Y
do
    grep -v "NA" chrom$chrom/*qassoc | grep -v "CH" | sort -g -k10 >> summary.pcutoff.1e-6
    echo $chrom
done
