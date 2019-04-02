import pandas as pd
import numpy as np 
import pdb 

#read in all the datasets 
accel_continuous=pd.read_csv("accelerometery_aggregate_phenotypes.continuous.txt",header=0,sep='\t',index_col=["FID","IID"])
print("loaded accel_continuous") 
exercise_slopes=pd.read_csv("Exercise_slopes.txt",header=0,sep='\t',index_col=["FID","IID"])
print("loaded exercise_slopes") 
hr_fitness=pd.read_csv("HR_fitnes.txt",header=0,sep='\t',index_col=["FID","IID"])
print("loaded hr_fitness") 
maxwd=pd.read_csv("MaxWD.txt",header=0,sep='\t',index_col=["FID","IID"])
print("loaded maxwd") 
recovery=pd.read_csv("Recovery.txt",header=0,sep='\t',index_col=["FID","IID"])
print("loaded recovery") 
accel_categorical=pd.read_csv("accelerometery_aggregate_phenotypes.categorical.txt",header=0,sep='\t',index_col=["FID","IID"])
print("loaded accel categorical") 

print("loaded phenotype dataframes") 

#merge continuous values 
continuous=accel_continuous.merge(exercise_slopes,left_index=True,right_index=True,how="outer")
continuous=continuous.merge(hr_fitness,left_index=True,right_index=True,how="outer") 
continuous=continuous.merge(maxwd,left_index=True,right_index=True,how="outer")
continuous=continuous.merge(recovery,left_index=True,right_index=True,how="outer") 
print("merged continuous values") 

#for continuous values, remove any outliers +/- 3std dev away from mean 
#temporarily set missing values to "NA" to allow for mean and std dev calculation 
continuous=continuous.astype('float32') 
continuous=continuous.replace(-1000,np.nan) 
continuous=continuous.replace(-9,np.nan)
print("replaced -1000 and -9 with NA") 
column_means=continuous._get_numeric_data().mean(axis=0,skipna=True)
print("calculated mean along axis 0") 
column_std=continuous._get_numeric_data().std(axis=0,skipna=True)
print("calculated std along axis 0")
continuous[continuous>(column_means+3*column_std)]=np.nan 
continuous[continuous<(column_means-3*column_std)]=np.nan 
print("removed outliers in continuous datasets > 3 std deviations away from the mean") 
print(str(continuous.shape)) 

#merge the datasets by subject id 
merged=continuous.merge(accel_categorical,left_index=True,right_index=True,how="outer")
print("added in categorical datasets") 

#replace missing values by -9 
merged=merged.replace(np.nan,-9) 
print("replaced na values w/ -9") 

#map the 22282 subject id's to the 24983 subject id's 
idmap=pd.read_csv("ukb22282_24983_mapping.tsv",header=None,sep='\t') 
idmap_dict=dict() 
for index,row in idmap.iterrows(): 
    id_22282=row[0] 
    id_24983=row[1] 
    idmap_dict[id_22282]=id_24983 
print("generated dictionary of id maps") 

#write updated phenotype file to disk  
outf=open("ashleylab_phenotypes_activity_and_fitness.txt",'w')
header=merged.columns 
outf.write("FID\tIID\t"+'\t'.join([str(i) for i in header])+'\n')
for index,row in merged.iterrows(): 
  cur_fid=index[0] 
  cur_iid=index[1] 
  try:
      translated_iid=idmap_dict[cur_iid] 
      translated_fid=translated_iid 
  except: 
      continue
  outf.write(str(translated_fid)+'\t'+str(translated_iid)+'\t'+'\t'.join([str(i) for i in row])+'\n')
