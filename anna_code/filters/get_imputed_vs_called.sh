#python get_imputed_vs_called.py --hits fitness_top_snps.csv \
#       --chrom_column 0 \
#       --pos_column 2 \
#       --snp_column 1 \
#       --calls_dir /oak/stanford/groups/euan/projects/ukbb/data/genetic_data/v2/EGAD00010001226/001 \
#       --outf fitness_top_snps_with_status.tsv

python get_imputed_vs_called.py --hits ../postprocessing/top_snps_activity_merged_plus_status.txt \
       --chrom_column 0 \
       --pos_column 2 \
       --snp_column 1 \
       --calls_dir /oak/stanford/groups/euan/projects/ukbb/data/genetic_data/v2/EGAD00010001226/001 \
       --outf ../postprocessing/top_snps_activity_merged_plus_status_updated.txt


