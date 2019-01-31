sbatch -J "merge_training" -o logs/merge.train.o -e logs/merge.train.e -p euan,owners  --time=48:00:00 --mem=100G  merge_plink_files_across_chroms.sh

