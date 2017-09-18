#for i in OverallAccelerationAverage TimeSpentOutdoorsSummer DurationWalkingForPleasure
for i in OverallAccelerationAverage
do
    for chrom in `seq 1 22`
    do
	sbatch -J "$i.$chrom.clump" -o logs/$i.$chrom.o -e logs/$i.$chrom.e -p euan,owners  clump.sh $i $chrom
    done
done
