for start_index in `seq 0 1000 96000`
do
    end_index=$(( $start_index + 1000))
    sbatch -J "transition$start_index" -o logs25/transition.$start_index.o -e logs25/transition.$start_index.e -p euan,owners get_number_transition_states.sh $start_index $end_index
done

    
