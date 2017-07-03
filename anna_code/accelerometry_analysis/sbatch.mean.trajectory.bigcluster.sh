export cluster_index=3
#for fragment in `seq 0 100 40461`
#for fragment in `seq 0 1000 23364`
for fragment in `seq 0 1000 18901`
do
    sbatch -J "mean_trajectory$cluster_index.$fragment" -o logs/mean_trajectory.$cluster_index.$fragment.o -e logs/mean_trajectory.$cluster_index.$fragment.e -p euan,owners extract_cluster_trajectories.sh $cluster_index $fragment
done
