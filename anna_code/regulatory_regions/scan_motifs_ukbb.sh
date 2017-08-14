#python scan_motifs_ukbb.py --pwm_list pwm_list_cisbp.txt --positions_bed activity_hits.30.bed --reference hg19/hg19.genome.fa --out_prefix activity_ --chrom_sizes hg19/hg19.chrom.sizes --num_hits_per_motif 1 --p_val 0.000001 --freqs --background_freqs activity_hits.foreground_freqs.txt --thresholds activity_hits.score_cutoffs.tsv  --totext

python scan_motifs_ukbb.py --pwm_list pwm_list_cisbp.txt --positions_bed fitness_hits.30.bed --reference hg19/hg19.genome.fa --out_prefix fitness_ --chrom_sizes hg19/hg19.chrom.sizes --num_hits_per_motif 1 --p_val 0.000001 --freqs --background_freqs fitness_hits.foreground_freqs.txt --thresholds fitness_hits.score_cutoffs.tsv  --totext

