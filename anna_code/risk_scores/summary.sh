# Subset the UKBB cohort to individuals with afib
generate_test_and_train_split.R
for chrom in `seq 1 22` 
do
    sbatch -J "$chrom.test" -o logs/$chrom.test.o -e logs/$chrom.test.e -p euan,owners  --time=48:00:00 --mem=20000  subset_genetics_to_afib.sh $chrom
done
sbatch.merge.sh


#Remove first degree relatives 
plink -bfile afib.cohort.train.1 --bmerge afib.cohort.test.1.bed afib.cohort.test.1.bim afib.cohort.test.1.fam  --make-bed --out afib.cohort.chrom1
sbatch -J kinship -o logs/kinship.o -e logs/kinship.e -p euan,owners -N 1 -c 16 --mem=20000 get_kinship.sh
python Filter_kinship.py
remove_first_degree_relatives.sh

# get SNP chromosome, position info from UKBB 
python get_chrom_pos.py
