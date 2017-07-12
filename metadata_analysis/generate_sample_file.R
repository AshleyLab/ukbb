#This document helps create the necessary sample file to work with BGEN format.  This can be used in snptest or plink1.9/2.
#In addition, the end of this document has a small portion that extracts the scores to create a "phenotype" file for plink1.9.


#Helper function: get_regex_cols
#Return the names of the columns of a particular matrix that match the desired pattern.

get_regex_cols<-function(cols,reg,...){return (cols[grepl(cols,pattern=reg,...)])}

#Note: Would be good to manually check if there is a discrepancy between calls and imputed again.

#Start by retrieving the IDs that are presented in the sample file for the BGEN files.
#We need to figure out what IDs have scores and imputed data.

complete_sample = read.table(file = "impv1.sample", header = T)
fam_sample = read.table(file = "chrom10.fam", header = F)
complete_sample = complete_sample[-1,]   #Remove the first line that contains the coding for the Oxford sample format.
complete_IDs = complete_sample[,1]
table(fam_sample$V1 %in% complete_IDs) # Check how many people in calls that are not in imputed (480).
table(complete_IDs %in% fam_sample$V1) # Check how many people in imputed are not in calls (2).
#fam_gen_exclusion = setdiff(fam_sample$V1, complete_IDs)
#fam_gen = setdiff(complete_IDs, fam_sample$V1)
#fam_gen_exclusion = c(fam_gen_exclusion, fam_gen) # The IDs that are not in either
relevant_IDs = as.integer(rownames(gwas_data)) # Get the IDs that we have scores for.
relevant_logical = relevant_IDs %in% complete_IDs  #Create a logical vector to extract the individuals that are scored and have imputed data.
excluded_relevant = relevant_IDs[!relevant_logical] # List of IDs that have scores but no imputed data.
relevant_IDs = relevant_IDs[relevant_logical] #List of IDs that have scores and imputed data.
table(excluded_relevant %in% fam_gen_exclusion) #the samples from the calls are not found in gen, all true
irrelevant_IDs = setdiff(complete_IDs, relevant_IDs) # Get the people that are genotyped but we don't have scores for.

#Retrieving the relevant information for the relevant ids.
#The first part generates corrected residuals with the covariates.

corrected_residuals = c()
is_residual = which(grepl("residuals",colnames(gwas_data),ignore.case = T))   # Get column names of residuals.
is_pcs = which(grepl("PC",colnames(gwas_data),ignore.case = T))   # Get column names of PCs.
for(ind in is_residual){
  y = gwas_data[,ind]
  x = gwas_data[,is_pcs]
  lm_obj = lm(y~x,na.action = na.exclude)
  new_residuals = lm_obj$residuals
  new_residuals = new_residuals[rownames(gwas_data)]
  names(new_residuals) = rownames(gwas_data)
  corrected_residuals = cbind(corrected_residuals,new_residuals)
}
rownames(corrected_residuals) = rownames(gwas_data)
colnames(corrected_residuals) = colnames(gwas_data)[is_residual]

#Start by creating the part of the sample file for the individuals with scores and imputed data.

missingness = get_regex_cols(colnames(pheno_data), "Missingness") # Get the column name of the missingness.
missingness = pheno_data[as.character(relevant_IDs), missingness] # Get the missingness.
sex = get_regex_cols(colnames(pheno_data), "Sex") # Get the column name of the sex.
sex = pheno_data[as.character(relevant_IDs), sex] # Get the sex.
gwas_sample = matrix(c(relevant_IDs,relevant_IDs, missingness, sex), nrow = length(relevant_IDs)) # Create a matrix following Oxford sample format: ID_1, ID_2, missingness.
rownames(gwas_sample) = as.character(relevant_IDs)
colnames(gwas_sample) = c("ID_1", "ID_2", "missing", "sex")
relevant_gwas_data = gwas_data[as.character(relevant_IDs),]  # Extract the scores provided by DD.
gwas_sample = cbind(gwas_sample, relevant_gwas_data) # Add the scores to our matrix.

#Follow with creating "NA" data for the rest of the individuals.

missingness = rep(NA, length(irrelevant_IDs)) # NAs for missingness column.
sex = get_regex_cols(colnames(pheno_data), "Sex") # Get the column name of the sex.
sex = pheno_data[as.character(irrelevant_IDs), sex] # Get the sex.
irrelevant_gwas_sample = matrix(c(irrelevant_IDs, irrelevant_IDs, missingness, sex), nrow = length(irrelevant_IDs))  # Create the second portion of the sample file.
colnames(irrelevant_gwas_sample) = c("ID_1", "ID_2", "missing", "sex")
dummy_vector = rep(NA, length(irrelevant_IDs)) #Again, a NA vector that has length of irrelevant_IDs
irrelevant_gwas_information = matrix(c(dummy_vector), nrow = length(irrelevant_IDs)) #Create a matrix with the NA vectors.
colnames(irrelevant_gwas_information) = colnames(relevant_gwas_data) 
irrelevant_gwas_sample = cbind(irrelevant_gwas_sample, irrelevant_gwas_information) # Append column wise to create the finished second portion of the sample file for no-score individuals.
rownames(irrelevant_gwas_sample) = as.character(irrelevant_IDs)

#We need to combine the two portions of the sample for the relevant and irrelevant IDs to create our finished sample file.

complete_gwas_sample = rbind(gwas_sample, irrelevant_gwas_sample) #Bind by row.
complete_gwas_sample = complete_gwas_sample[as.character(complete_IDs),, drop = FALSE] #Reorder to the original order!!!
write.table(irrelevant_IDs, quote = F, col.names = F, row.names = F, file = "exclusionmay17.txt")  #List of IDs for exclusion (no score but have imputed data).
#The previous line is necessary for snptest, but is most likely unnecessary for plink.  In plink1.9, use the option --prune to remove the individuals without
#phenotypes.  For plink2, use --require-pheno.
val = c(0,0,0) #Coding information for the Oxford sample format.
val1 = c("D","P") #Coding information for the Oxford sample format.
#val2 = c("D","D","C","D","C","D","D","D","D","C","C","D","D","D","D","D","D","D","D","D","D","D","D","D","C","D","C","D","C","C","C")
val = c(val, val1) #Combine the vectors.
filename = "testjune1.sample" #Name the sample file.
write(gsub(" ", "_",colnames(complete_gwas_sample)), file = filename, ncol = length(colnames(complete_gwas_sample))) #Write the column names first as a header.
write(val, file = filename, append = TRUE, ncol = length(colnames(complete_gwas_sample))) #Then append the mandatory Oxford coding.
write.table(complete_gwas_sample, file = filename, append = TRUE, col.names = F, row.names = F, quote = F) #Finish by adding the rest of the data in the matrix.

#This portion makes a file of scores with IDs for the plink program, which is to be used with --pheno for plink1.9.

fname = "x" 
residual_columns = get_regex_cols(colnames(complete_gwas_sample), "Residuals")
plink_header = c("FID", "IID" ,gsub(" ", "_", residual_columns)
write(plink_header, file = fname, ncol = length(plink_header))
write.table(complete_gwas_sample[,c("ID_1", "ID_2", residual_columns)], file = fname, append = T, col.names = F, row.names = F, quote = F)