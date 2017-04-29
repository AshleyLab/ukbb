#!/bin/bash

map=~/projects1/ukbb/file_handlers/genotype_map.csv
out=~/projects1/ukbb/data/genetic_data/calls/plink/
filehandler=~/projects1/ukbb/file_handlers
in=~/projects1/ukbb/data/genetic_data/calls

# loop over files
for i in `ls  ~/projects1/ukbb/data/genetic_data/calls/chrom*`;
do
    ped=`basename $i .cal`
    cmd1="$filehandler/gconv $in/$i ${out}/${ped} ped -m${map}"
    cmd2="plink --file ${out}/${ped} --make-bed --out ${out}/${ped}"
    echo $out/$ped
    echo '#!/bin/bash' > run.sh
    echo $cmd2 >> run.sh
    sbatch -J "conv-$i" -e logs/${ped}.err -o logs/${ped}.out -p euan,owners run.sh
    #sbatch -J "conv-$i" -p euan,owners run.sh
    #exit
done

