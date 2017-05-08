kmeansAIC = function(fit){
  
  m = ncol(fit$centers)
  n = length(fit$cluster)
  k = nrow(fit$centers)
  D = fit$tot.withinss
  return(D + 2*m*k)
}

kmeansBIC = function(fit){
  
  m = ncol(fit$centers)
  n = length(fit$cluster)
  k = nrow(fit$centers)
  D = fit$tot.withinss
  return(D + log(n)*m*k)
}


#cdf 
cdf<-function(x)
{
  x<-x/sum(x)
  cdfval=c()  
  total_sum=0 
  for(i in 1:length(x))
  {
    total_sum=total_sum+x[i] 
    cdfval=c(cdfval,total_sum)
  }
  return(cdfval)
}

#sum of squares 
ss<-function(x)
{
  sum(x^2)
}

#normalization 
normalize<-function(x)
{
  (x-min(x))/(max(x)-min(x))
}

#takes a difference between the last and first timestamp and divides by the number of samples to estimate the sampleing frequency
get_fs<-function(x)
{
  if(typeof(x)=="double")
  {
    delta<-x[length(x)]-x[1]
  }
  else
  {
    #convert from string to timestamp
    last <- strptime(x[length(x)],"%Y-%m-%dT%H:%M:%S")
    first<-strptime(x[1],"%Y-%m-%dT%H:%M:%S")
    delta<-as.numeric(difftime(last,first,units="secs"))
  }
  fs=length(x)/delta
  return(c(fs=fs,delta=delta)) 
}

#MMAV--modified mean absolute value 
get_mmav<-function(x)
{
  N=length(x)
  mmav=0 
  for(i in 1:N)
  {
    if(i<0.25*N)
    {
      Wn=0.5
    }
    else if(i>0.75*N)
    {
      Wn=0.05
    }
    else
    {
      Wn=1
    }
    mmav=mmav+Wn*abs(x[i])
  }
  mmav=mmav/N
  return(mmav)
}

get_zero_crossing<-function(x)
{
  lowgroup=x[1:length(x)-1]
  highgroup=x[2:length(x)]
  product=sign(lowgroup*highgroup)
  positive=which(product %in% 1)
  product[positive]=0 
  delta=sign(lowgroup-highgroup)
  negative=which(delta %in% -1)
  delta[negative]=0 
  zero_val=abs(sum(product*delta))
  return(zero_val)
}

#frequency ratio: A measure termed the frequency
#ratio is defined as the ratio of the power in the high frequency band from 8 to 24 Hz
#compared to the power in the low frequency band from 3 to 5 Hz
get_frequency_ratio<-function(s,low_band_low,low_band_high,high_band_low,high_band_high)
{
  low_band_low_index<-max(which(s$freq<low_band_low))+1 
  low_band_high_index<-max(which(s$freq<low_band_high))+1 
  high_band_low_index<-max(which(s$freq<high_band_low))+1 
  high_band_high_index<-max(which(s$freq<high_band_high))+1 
  low_freq_spectrum<-sum(s$spec[low_band_low_index:low_band_high_index])
  high_freq_spectrum<-sum(s$spec[high_band_low_index:high_band_high_index])
  ratio<-low_freq_spectrum/high_freq_spectrum 
  return(ratio)
}

#calculates the period of a signal 
get_period<-function(x,interval)
{ 
  x_sub<-x[round(length(x)/2-interval/2):round(length(x)/2+interval/2)]
  x_unique=unique(x_sub) 
  if (length(x_unique)<0.5*length(x_sub))
  {
  #We have excessive duplicates, remove them 
  x_sub=x_unique 
  }	
  #is most of the signal power above the mean (don't flip for peak calls) or below the mean (flip for peak calls)
  x_mean<-mean(x_sub)
  x_high<-subset(x_sub,x_sub>x_mean)
  x_high_mean<-mean(x_high)
  x_low<-subset(x_sub,x_sub<x_mean)
  x_low_mean<-mean(x_low)
  if (abs(x_high_mean)>abs(x_low_mean))
  {
    peaks<-findPeaks(x_sub)
    period<-mean(diff(peaks))
    for (i in 1:(length(peaks)-1))
    {
    delta<-abs(peaks[i+1]-peaks[i])
    if (min(delta,period)/max(delta,period)>0.5)
    {
       return(x[peaks[i]:peaks[i+1]])
    }
    }
  }
  else
  {
    peaks<-findPeaks(-1*x_sub) 
    period<-mean(diff(peaks))
    for(i in 1:(length(peaks)-1))
    {
      delta<-abs(peaks[i+1]-peaks[i])
      if(min(delta,period)/max(delta,period)>0.5)
      {
        return(x[peaks[i]:peaks[i+1]])
      }
    }
    
  }
  return(NULL)
}
