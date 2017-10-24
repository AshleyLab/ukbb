sbatch -J "aggregate_clumped_for_23andme" -o logs/aggregate_for_23andme.o -e logs/aggregate_for_23andme.e -p akundaje,owners --mem=20000 aggregate_clump.sh $task 
