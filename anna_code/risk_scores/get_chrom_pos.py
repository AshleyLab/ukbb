import pandas as pd 
as_data="MEGASTROKE_data/MEGASTROKE.4.CES.EUR.out"
bim_data="afib.cohort.train.bim"

as_data = pd.read_csv(as_data, header = 0, sep= " ")
bim = pd.read_csv(bim_data, header = None, sep= "\t")
bim.columns=['chr','rs','info','pos','a1','a2' ]
as_data.columns=['rs','alt','ref','altfreq','efalt','std','pval']
merged=as_data.merge(bim,left_on='rs',right_on='rs')
merged.to_csv("Merged.CES.EUR.tsv",sep='\t')


