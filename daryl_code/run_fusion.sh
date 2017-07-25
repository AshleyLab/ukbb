#!/bin/bash
module load R/3.3.3.gcc
gwass=`ls /oak/stanford/groups/euan/projects/ukbb/da_dh/gwas/immediate_gwas2/*linear`
#gwas=/oak/stanford/groups/euan/projects/ukbb/da_dh/gwas/jul17_hr_pred_100_plink2/gwas2/chr16impv1.Residuals_HR_pred_100_plink.glm.linear
phenoname=Residuals_HR_pred_100
tissues=`ls /share/PI/euan/apps/fusion_twas-master/WEIGHTS/*pos`
#tissues=/share/PI/euan/apps/fusion_twas-master/WEIGHTS/GTEx.Muscle_Skeletal.pos
out=/oak/stanford/groups/euan/projects/ukbb/twas/

# munge 
echo $chr

# loop over phenotypes/gwas
# loop over chr
# loop over weights/tissues
for gwas in $gwass; do
    gwasname=`basename $gwas`
    chr=`echo $gwasname | sed -e 's/^chr//' -e 's/impv1.Residuals_HR_pred_100_plink.glm.linear//'`
    echo "SNP A1 A2 N Z" > $out/tmp/${gwasname}.tmp 
    #awk 'NR>1{print $3,$4,$5,$7,$10}' $gwas >> $out/tmp/${gwasname}.tmp 
    awk 'NR>1{print $2,$4,$8}' $gwas >> $out/tmp/${gwasname}.tmp 
    gwas=$out/tmp/${gwasname}.tmp
    
    for tissue in $tissues; do
        tissuename=`basename $tissue | sed -e 's/GTEx.//' -e 's/.pos//'`
        Rscript ~/apps/fusion_twas-master/FUSION.assoc_test.R \
            --sumstats $gwas   \
            --weights $tissue \
            --weights_dir /share/PI/euan/apps/fusion_twas-master/WEIGHTS \
            --ref_ld_chr /share/PI/euan/apps/fusion_twas-master/LDREF/1000G.EUR. \
            --chr $chr \
            --out $out/${gwasname}_${tissuename}.dat
    done
done
