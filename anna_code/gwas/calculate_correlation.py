import pandas as pd
data =pd.DataFrame.from_csv("accelerometery_aggregate.txt",header=0,sep='\t')
pearson_corval= data.corr(method='pearson', min_periods=1)
pearson_corval=pearson_corval.round(3)
pearson_corval.to_csv(path="pearson_correlation.tsv",sep='\t')
spearman_corval= data.corr(method='spearman', min_periods=1)
spearman_corval=spearman_corval.round(3)
spearman_corval.to_csv(path="spearman_correlation.tsv",sep='\t')

