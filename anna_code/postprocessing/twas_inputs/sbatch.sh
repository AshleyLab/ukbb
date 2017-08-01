for i in `seq 19 22`
do
    sbatch -J "make_twas_input_$i" -o logs/make_twas_input_$i.o -e logs/make_twas_input_$i.e -p euan,owners --mem=15000 make_twas_inputs.sh $i
done

