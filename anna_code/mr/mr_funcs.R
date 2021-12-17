rm(list=ls())
library('MendelianRandomization')
args <- commandArgs(TRUE)
#exposure_name=args[1]
#outcome_name=args[2]
exposure_name="OverallAccelerationAverage"
#exposure_name="X12am_6am"
outcome_name="CvdStatus"

covs=read.table("mr_covariates.txt",header=TRUE,row.names=1,sep='\t')
covs$IID=NULL
covs$Sex=factor(covs$Sex)
#covs$f.batch=factor(covs$f.batch)
covs$f.batch=NULL

snps=read.table(paste("mr_snp_subsets/",exposure_name,sep=""),header=TRUE,sep='\t',row.names=1)
snps[snps==-1]=NA

exposure_mr=read.table("mr_exposures.txt",header=TRUE,row.names=1,sep='\t')
expo_v=as.matrix(subset(exposure_mr,select=c(exposure_name)))

outcome_mr=read.table("mr_outcome.txt",header=TRUE,row.names=1,sep='\t')
outcome_v=factor(as.matrix(subset(outcome_mr,select=c(outcome_name))))

#expo_v=expo_v[is.na(expo_v)==FALSE]
#outcome_v=outcome_v[is.na(expo_v)==FALSE]
#snps=snps[is.na(expo_v)==FALSE]
#covs=covs[is.na(expo_v)==FALSE,]

get_lm_stats <- function(x, y) {
  if (length(unique(y[!is.na(y)])) == 2) {
    o = glm(y ~ x, family = binomial(link = 'logit'))
    o = summary(o)$coefficients
  }
  else{
    o = summary(lm(y ~ x))$coefficients
  }
  b = o[-1, 1]
  s = o[-1, 2]
  p = o[-1, 4]
  return(c(b, s, p))
}

get_lm_stats_with_covs <- function(x, y, covs) {
  if (is.null(covs)) {
    return (get_lm_stats(x, y))
  }
  names(x) = c("x")
  names(y) = c("y")
  d = data.frame(x, y, covs)
  #check whether the y is a factor (logistic regression) or a continuous value (linear regression)
  if (is.factor(y) == TRUE)
  {
    o = summary(glm(y ~ ., data = d, family = binomial))$coefficients
  }
  else
  {
    o = summary(lm(y ~ ., data = d))$coefficients
  }
  b = o["x", 1]
  s = o["x", 2]
  p = o["x", 4]
  return(c(b, s, p))
}

min_effect_size = 0
    covs = covs[rownames(snps), ]
    #we make the assumption that expo_snps and outcome_snps are the same
    
    lm_res = apply(snps,
                   2,
                   get_lm_stats_with_covs,
                   y = expo_v,
                   covs = covs)
    bx = lm_res[1, ]
    bxse = lm_res[2, ]
    lm_res = apply(snps,
                   2,
                   get_lm_stats_with_covs,
                   y = outcome_v,
                   covs = covs)
    by = lm_res[1, ]
    byse = lm_res[2, ]
    print("Two sample estimate completed")
    #browser()     
    to_keep = abs(bx) >= min_effect_size
    if (sum(to_keep) == 0) {
      to_keep = abs(bx) >= median(abs(bx))
    }
    snps[is.na(snps)] = 0
    snps = as.matrix(snps)
    
    weighted_causal_v = as.matrix(snps[, to_keep]) %*% as.matrix(bx[to_keep])
    if (length(unique(expo_v[!is.na(expo_v)])) == 2) {
      lm_expo_vs_g = glm(expo_v ~ weighted_causal_y, family = binomial)
    }else{
      lm_expo_vs_g = lm(expo_v ~ weighted_causal_v)
    }
    bx_obj = summary(lm_expo_vs_g)$coefficients
    bx_u = bx_obj[-1, 1]
    bxse_u = bx_obj[-1, 2]
    
    if (length(unique(outcome_v[!is.na(outcome_v)])) == 2) {
      lm_outcome_vs_g = glm(outcome_v~ weighted_causal_v, family = binomial)
    }else
    {
      lm_outcome_vs_g = lm(outcome_v ~ weighted_causal_v)
    }
    by_obj = summary(lm_outcome_vs_g)$coefficients
    by_u = by_obj[-1, 1]
    byse_u = by_obj[-1, 2]
    
    res = list()
    if (sum(to_keep) > 1) {
      mr_in = mr_input(bx[to_keep], 
                       bxse[to_keep], 
                       by[to_keep], 
                       byse[to_keep],
                       snps=colnames(snps),
                       exposure=exposure_name,
                       outcome=outcome_name)
      #res[["multivar_input"]] = mr_in
      res[["multivar_all"]] = mr_allmethods(mr_in)
      #MRAllObject_all <- mr_allmethods(mr_in, method = "all")
      png(paste(exposure_name,outcome_name,'png',sep='.'),height=10,width=10,units='in',res=600)
      mr_plot(mr_in,interactive=FALSE,labels=TRUE)
      dev.off() 
      mr_egger(mr_in, F, T)
      mr_median(mr_in)
    }else{
      mr_in = mr_input(bx_u, 
                       bxse_u, 
                       by_u, 
                       byse_u,
                       snps=colnames(snps),
                       exposure=exposure_name,
                       outcome=outcome_name)
      #res[["univar_input"]] = mr_in
      #res[["univar_ivw"]] = mr_ivw(mr_in)
      #res[["univar_ml"]] = mr_maxlik(mr_in)
      #res[["univar_expo_lm"]] = lm_expo_vs_g
      #res[["unival_outcome_lm"]] = lm_outcome_vs_g
      #MRAllObject_all <- mr_allmethods(mr_in, method = "all")
      #mr_plot(MRAllObject_all)
      png(paste(exposure_name,outcome_name,'png',sep='.'),height=10,width=10,units='in',res=600)
      mr_plot(mr_in,interactive=FALSE,labels=TRUE)
      dev.off() 
    }
    