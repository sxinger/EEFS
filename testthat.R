## testing script ##
rm(list=ls()); gc()

setwd("~/#R-pkg/EEFS")

library("dplyr")
library("tidyr")
library("magrittr")

#===test feature_ensemble() function
source("./R/feature_ensemble.R")
#testing parameters
feat_rk1<-data.frame(Feature=c("v1","v2","v3","v4","v5"),
                     Rank=1:5)

feat_rk2<-data.frame(Feature=c("v1","v2","v3","v4"),
                     Rank=c(4,1,2,3))

#test 1: without weighting
rank_ens<-feature_ensemble(
  rank_lst=list(feat_rk1,feat_rk2),
  var_colnm="Feature",
  rank_colnm="Rank",
  ensemble_mth<-"exp_k",
  k=3
)


rank_ens<-feature_ensemble(
  rank_lst=list(feat_rk1,feat_rk2),
  var_colnm="Feature",
  rank_colnm="Rank",
  ensemble_mth="exp_k_wt",
  k=3
)
