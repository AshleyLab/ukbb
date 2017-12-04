for f in `seq 0 784`
do
    sbatch -J "average_intervas$f" -o logs/average_intervals.$f.o -e logs/average_intervals$f.e -p euan,owners get_interval_averages.sh $f
done
