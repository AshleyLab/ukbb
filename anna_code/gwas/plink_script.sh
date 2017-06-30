export chrom=chrom$1
plink --bfile /scratch/PI/euan/projects/ukbb/data/genetic_data/calls/plink/$chrom --pheno accelerometry.pheno --assoc  --all-pheno --pfilter 1e-6 --missing-phenotype -1000 --out $chrom



