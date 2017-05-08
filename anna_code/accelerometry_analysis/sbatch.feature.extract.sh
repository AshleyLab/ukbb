#958
for subject_set in `seq 897 958`
do
    sbatch -J "feat.extract.$subject_set" -o logs/feat.extract.$subject_set.o -e logs/feat.extract.$subject_set.e -p euan,owners feat.extract.sh $subject_set 
done
