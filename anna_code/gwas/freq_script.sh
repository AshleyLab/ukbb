for chrom in `seq 1 22` X Y
do
    plink --bfile /scratch/PI/euan/projects/ukbb/data/genetic_data/calls/plink/chrom$chrom --freq --out /scratch/PI/euan/projects/ukbb/data/genetic_data/freqs/chrom$chrom
done
