#958
for i in `seq 1 784`
do
    sbatch -J "feat.extract.$i" -o logs/feat.extract.$i.o -e logs/feat.extract.$i.e -p euan,owners extract_features_small.sh $i
done
