## testing script ##
rm(list=ls()); gc()

setwd("~/#R-pkg/EEFS")

library("dplyr")
library("tidyr")
library("magrittr")
library("pROC")
library("ROCR")

#===simulate data
feat_rk1<-data.frame(Feature=c("v1","v2","v3","v4","v5"),
                     Rank=1:5)

feat_rk2<-data.frame(Feature=c("v1","v2","v3","v4"),
                     Rank=c(4,1,2,3))

feat_rk3<-data.frame(Feature=c("v1","v6","v7","v8"),
                     Rank=c(4,1,2,3))

#===test feature_ensemble() function
source("./R/feature_ensemble.R")
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
  wt=c(0.8,0.7),
  k=3
)


#===test evaluation module
source("./R/evaluation.R")
#test 1: eval_stability()
eval_stability(rank_lst=list(feat_rk1,feat_rk2),
               var_colnm="Feature",
               rank_colnm="Rank",
               metric="kuncheva",
               f=10,
               d=5)


eval_stability(rank_lst=list(feat_rk1,feat_rk3),
               var_colnm="Feature",
               rank_colnm="Rank",
               metric="kuncheva",
               f=10,
               d=5)


eval_stability(rank_lst=list(feat_rk1,feat_rk2),
               var_colnm="Feature",
               rank_colnm="Rank",
               metric="wci",
               f=10,
               d=5)


eval_stability(rank_lst=list(feat_rk1,feat_rk3),
               var_colnm="Feature",
               rank_colnm="Rank",
               metric="wci",
               f=10,
               d=5)



eval_stability(rank_lst=list(feat_rk1,feat_rk2),
               var_colnm="Feature",
               rank_colnm="Rank",
               metric="wci_rel",
               f=10,
               d=5)


eval_stability(rank_lst=list(feat_rk1,feat_rk3),
               var_colnm="Feature",
               rank_colnm="Rank",
               metric="wci_rel",
               f=10,
               d=5)

