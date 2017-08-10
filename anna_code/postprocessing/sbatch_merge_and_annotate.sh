sbatch -J "merge_and_annotate" -o logs/merge_and_annotate.o -e logs/merge_and_annotate.e -p euan,owners --mem=20000 merge_and_annotate.sh

