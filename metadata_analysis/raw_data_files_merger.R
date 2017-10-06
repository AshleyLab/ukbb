try({setwd("/Users/davidhsu/Documents/ukbb")})
try({setwd("/Users/david/Desktop/ukbb/")})
try({setwd('/scratch/PI/euan/projects/ukbb/da_dh/')})
try({setwd('/scratch/PI/euan/projects/ukbb/david/')})
# We cannot change the working directory as we work in RStudio, so we comment out this line temporarily.
trsource("auxiliary_functions.R")

# This script collates all pheno data files
# This part both merges the different initial pheno data files that
# we got from the uk biobank scripts (e.g., ukb7454.r) and
# translates the column names into human readable names
# we save these preprocessing results in a new data frame - "biobank_collated_pheno_data.RData"
# so you can just load it and skip  ahead
# We also exclude subjects that did not do the clinical test
load('biobank_pheno_data.RData')
load('biobank_pheno_data_addition_1.RData')
pheno_data = as.data.frame(pheno_data)
table(is.na(pheno_data))
pheno_data_addition_1 = as.data.frame(pheno_data_addition_1)
# Step 1: put it all in one big data frame
# Some tests - all should be TRUE
dim(pheno_data);dim(pheno_data_addition_1)
# all(pheno_data[,1]==pheno_data_addition_1[,1]) # Is the subject order the same? Answer: no!
length(setdiff(pheno_data[,1],pheno_data_addition_1[,1]))==0
length(unique(pheno_data[,1])) == nrow(pheno_data)
rownames(pheno_data) = pheno_data[,1]
rownames(pheno_data_addition_1) = pheno_data_addition_1[,1]
pheno_data_addition_1 = pheno_data_addition_1[rownames(pheno_data),]
all(rownames(pheno_data)==rownames(pheno_data_addition_1)) # Is the subject order the same? Answer: now yes
# Are all the columns in the addition new? only the id should appear in the intersection
intersect(colnames(pheno_data_addition_1),colnames(pheno_data)) == 'f.eid'
# merge the tables
pheno_data = cbind(pheno_data,pheno_data_addition_1[,-1])
dim(pheno_data)
rm(pheno_data_addition_1)
# remove the first column - we have this information as the rownames
pheno_data = pheno_data[,-1]

# Load the field names
# We want informative column names - makes it easier
# for analysis with regular expressions
field_id2name_raw = read.delim('field_id_to_name.txt',sep='\t')
field_id2name = as.character(field_id2name_raw[,2])
names(field_id2name) = as.character(field_id2name_raw[,1])
raw_colnames = colnames(pheno_data)
raw_colnames_table = t(sapply(raw_colnames,function(x) strsplit(x,split='\\.+')[[1]]))
# the second column in the matrix that we created is the field id
# transform the id to the names using the vector above
raw_colnames_table[,2] = field_id2name[raw_colnames_table[,2]]
new_col_names = apply(raw_colnames_table[,-1],1,paste,collapse = '.')
colnames(pheno_data) = new_col_names
# exclude samples without time series of the clinical test
# our criteria is to check if all workload data are NAs then
# the subject should be excluded
time_regex = "ECG, phase time"
workload_regex = "ECG, load"
time_columns = colnames(pheno_data)[grepl(colnames(pheno_data),pattern = time_regex)]
workload_columns = colnames(pheno_data)[grepl(colnames(pheno_data),pattern = workload_regex)]
time_col_mat = pheno_data[,time_columns]
time_points_is_all_na_in_row = apply(is.na(time_col_mat),1,all)
table(time_points_is_all_na_in_row)
wld_col_mat = pheno_data[,workload_columns]
wld_points_is_all_na_in_row = apply(is.na(wld_col_mat),1,all)
# sanity check: whenever we have a workload score we should have a time point
length(setdiff(which(!is.na(wld_col_mat)),which(!is.na(time_col_mat))))==0
# Sanity check: the subject sets should be similar
table(wld_points_is_all_na_in_row,time_points_is_all_na_in_row) #4/13/2017: the sets agree except for 196 subjects with time points and no worlkloads
rowSums(!is.na(time_col_mat[which(!time_points_is_all_na_in_row & wld_points_is_all_na_in_row)[1:10],]))
subject_cleaned_pheno_data = pheno_data[!wld_points_is_all_na_in_row,]
save(pheno_data,file="biobank_collated_pheno_data.RData")
save(subject_cleaned_pheno_data,file="biobank_collated_filtered_pheno_data.RData")
