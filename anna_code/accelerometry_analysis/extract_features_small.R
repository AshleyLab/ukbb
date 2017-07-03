library(wavelets) 
dwt_transform_features_small <- function(data,nlevels,a,b)
{
d_energy_overall=0 
d_energy_by_level=c()
smv=sum(data*data)
dwt_x=dwt(data,filter='d10',n.level=nlevels)
  for (i in 1:nlevels)
  {
    #GET THE TOTAL ENERGY AT THE LEVEL
    d_energy = 0
    for (j in 1:i)
    {
      d_energy_component = sum((dwt_x@W[[j]]) *(dwt_x@W[[j]]))
      d_energy = d_energy + d_energy_component
    }
  d_energy_overall=d_energy_overall+d_energy 
  d_energy_by_level=c(d_energy_by_level,d_energy)
  }
feat1=(d_energy_by_level[a]+d_energy_by_level[b])/d_energy_overall 
feat2=(d_energy_by_level[a]+d_energy_by_level[b])/smv 
return(c(feat1,feat2))
}
#datafile="/home/anna/ukbb/dwt/examples/1129255_90004_0_0.csv"
args <- commandArgs(TRUE)
datafile=args[1]
outfile=args[2]
#print(datafile)
data=as.matrix(read.table(datafile,header=FALSE,sep=','))
parts=strsplit(datafile,'/')[[1]]
subject=strsplit(parts[length(parts)],'_')[[1]][1]
result_dwt<-dwt_transform_features_small(data,8,5,6)
out_string=c(subject,result_dwt[1],result_dwt[2])
write.table(t(out_string),file=args[2],quote=FALSE,append=TRUE,sep = '\t',row.names=FALSE,col.names=FALSE)
