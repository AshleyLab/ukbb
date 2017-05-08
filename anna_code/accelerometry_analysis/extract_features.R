#EXTRACTS FEATURES FROM A TIMESERIES SIGNAL, RETURNS A DATAFRAME

#RESOURCES 
#FFT & TIME SERIES 
#http://stats.stackexchange.com/questions/50807/features-for-time-series-classification
#https://github.com/nickgillian/grt/blob/master/GRT/FeatureExtractionModules/FFT/FFT.h
#https://cran.r-project.org/web/views/TimeSeries.html
library(wavelets) 

##############################################################################################################
dwt_transform_features <- function(x,subject,nlevels)
{
  #TAKE THE DWT TRANSFORM
  #Get the filter and the number of levels from the Parameters file
  #To clarify "W" are the wavelet (detail) coefficients; "V" are the scaling (approximation) coefficients
  dwt_x=dwt(x,filter = 'd4',n.level = nlevels)
  #Extract a feature vector at each of the levels
  feature_vectors=list() 
  for (i in 1:nlevels)
  {
    #GET THE TOTAL ENERGY AT THE LEVEL
    a_energy = sum((dwt_x@V[[i]]) ^ 2)
    d_energy = 0
    d_energy_components = c()
    EDR_d_names=c() 
    for (j in 1:i)
    {
      d_energy_component = sum((dwt_x@W[[j]]) ^ 2)
      d_energy_components = c(d_energy_components,d_energy_component)
      d_energy = d_energy + d_energy_component
      EDR_d_names=c(EDR_d_names,paste("EDR_d",j,sep=""))
    }
    total_energy = a_energy + d_energy
    
    #GET THE ENERGY RATIOS AT THE LEVEL (Ayrulu-Erdem & Barshan, doi:10.3390/s110201721)
    EDR_a = a_energy / total_energy
    EDR_d = d_energy_components / total_energy
    
    #GET THE FEATURES INCORPORATING ENERGY AND COEFFICIENTS
    min_EDR_d = min(EDR_d)
    max_EDR_d = max(EDR_d)
    mean_EDR_d = mean(EDR_d)
    if(length(d_energy_components)==1)
    {
      var_EDR_d=0 
    }
    else
    {
      var_EDR_d = var(EDR_d)
    }
    #normalized means of decomposition coefficients
    means = c(mean(dwt_x@V[[i]]))
    mean_names=c("mean_scaled")
    for (j in 1:i)
    {
      means = c(means,mean(dwt_x@W[[j]]))
      mean_names=c(mean_names,paste("means_wav",j,sep=""))
    }
    #means = normalize(means)
    
    #variances of decomposition coefficients
    vars = c(var(dwt_x@V[[i]]))
    var_names=c("var_scaled")
    for (j in 1:i)
    {
      vars = c(vars,var(dwt_x@W[[j]]))
      var_names=c(var_names,paste("var_wav",j,sep=""))
    }
    #vars = normalize(vars)
    
    #combine into a dataframe 
    level_features<-data.frame(t(c(total_energy,EDR_a,EDR_d,min_EDR_d,max_EDR_d,mean_EDR_d,var_EDR_d,means,vars)),row.names=subject)
    names(level_features)<-(c("TotalEnery","EDR_a",EDR_d_names,"min_EDR_d","max_EDR_d","mean_EDR_d","var_EDR_d",mean_names,var_names))
    feature_vectors[[i]]=level_features
  }
  return(feature_vectors) 
}
