sbatch -J "aggregate_covariates1" -o aggregate1.o -e aggregate1.e -p euan,owners --time=04:00:00 --mem=15000 aggregate_covariates.sh
#sbatch -J "aggregate_fields" -o aggregate.o -e aggregate.e -p euan,owners --time=04:00:00 --mem=15000 aggregate_fields.sh
