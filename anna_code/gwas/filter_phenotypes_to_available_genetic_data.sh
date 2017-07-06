#Quantitative 
#python filter_phenotypes_to_available_genetic_data.py --input_df /scratch/PI/euan/projects/ukbb/gwas/accelerometry_plink/v1/accelerometery_aggregate_phenotypes.continuous.txt \
#       --subject_ids /scratch/PI/euan/projects/ukbb/data/genetic_data/subject_ids.txt \
#       --outf /scratch/PI/euan/projects/ukbb/gwas/accelerometry_plink/v1/accelerometery_aggregate_phenotypes.continuous.filtered.txt

#Categorical 
#python filter_phenotypes_to_available_genetic_data.py --input_df /scratch/PI/euan/projects/ukbb/gwas/accelerometry_plink/v1/accelerometery_aggregate_phenotypes.categorical.txt \
#       --subject_ids /scratch/PI/euan/projects/ukbb/data/genetic_data/subject_ids.txt \
 #      --outf /scratch/PI/euan/projects/ukbb/gwas/accelerometry_plink/v1/accelerometery_aggregate_phenotypes.categorical.filtered.txt


#Filter covariates
python filter_phenotypes_to_available_genetic_data.py --input_df /scratch/PI/euan/projects/ukbb/gwas/accelerometry_plink/v1/accelerometery_covariates.txt \
       --subject_ids /scratch/PI/euan/projects/ukbb/data/genetic_data/subject_ids.txt \
       --outf /scratch/PI/euan/projects/ukbb/gwas/accelerometry_plink/v1/accelerometery_covariates.filtered.txt
