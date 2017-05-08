for subject_set in `seq 0 1036`
do
    sbatch -J "timealign_$subject_set" -o logs/timealign.$subject_set.o -e logs/timealign.$subject_set.e -p euan,owners time_align_accelerometry_data.sh $subject_set 
done
