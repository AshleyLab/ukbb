sbatch -J "aggregate_clumped_for_23andme" -o logs/aggregate_for_23andme.o -e logs/aggregate_for_23andme.e -p euan,owners --mem=20000 summarize_clumped_data_for_23andme.sh 
