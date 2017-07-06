sbatch -J "aggregate_accel_pheno" -o aggregate.o -e aggregate.e -p euan,owners --time=04:00:00 --mem=15000 aggregate_fields.sh
