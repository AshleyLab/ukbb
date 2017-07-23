for i in `seq 2 233` 
do
    sbatch -J "$i.act" -o logs/$i.act.o -e logs/$i.act.e -p akundaje,owners --mem 5000  perform_icd_phenotype_regression.activity.sh $i 
done

for i in `seq 2 11`
do
sbatch -J "$i.act" -o logs/$i.act.o -e logs/$i.act.e -p akundaje,owners --mem 5000  perform_icd_phenotype_regression.other.sh $i 
done
for i in `seq 2 9`
do
    sbatch -J "$i.act" -o logs/$i.act.o -e logs/$i.act.e -p akundaje,owners --mem 5000  perform_icd_phenotype_regression.fitness.sh $i 
done 
