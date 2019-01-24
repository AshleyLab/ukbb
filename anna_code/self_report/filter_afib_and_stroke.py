#We want to filter file "atrial.fibrillation.1471.tsv" to individuals who have atrial fibrillation (1471) 
#or the combination of stroke (1081) but not brain hemmorhage (1491) or subarachnoid hemmorhage (1086) 
import pandas as pd 
import pdb 

#read in the data 
data=pd.read_csv("atrial.fibrillation.1471.tsv",header=0,sep='\t')

#extract individuals who had stroke but not brain hemmorhage or subarachnoid hemmorhage 
stroke_subset=data[data['1081']-data['1086']-data['1491']>0]

#extract individuals with afib 
afib_subset=data[data['1471']>0]
 
#merge the afib and stroke subset; drop the hemmorhage columns as we no longer need them 
afib_or_stroke=afib_subset.merge(stroke_subset,how="outer").drop(['1086','1491'],axis=1) 

#write output file 
afib_or_stroke.to_csv("afib_or_stroke_filtered_for_hemmorhage.tsv",sep='\t',index=False)

