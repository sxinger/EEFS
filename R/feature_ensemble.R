################################################################
# Feature Ensemble Module                                      #
# Parameters:                                                  #
# -- rank_lst:feature ranking list                             #
#    - list: each list can be a data.frame or matrix           #
#       with two columns (var_colnm, rank_colnm)               #
# -- var_colnm: column name for variable in each list          #
# -- rank_colnm: column name for rank in each list             #
# -- ensemble_mth: feature ranking ensemble methods            #
#    - mean: unweighted average of ranks                       #
#    - top_k: hard membership of a feature ranking among top k #                                 
#    - exp_k: soft membership of a  feature ranking among top k#                 
#    - wt_mean: mean adjusted by weights                       #
#    - wt_top_k: top_k adjusted by weights                     #
#    - wt_exp_k: exp_k adjusted by weights                     #
# -- wt: weights associated with each feature set '            #
#    - each row of rank_mt                                     #
# -- k: top k                                                  #
################################################################

feature_ensemble<-function(rank_lst,
                           var_colnm="var",
                           rank_colnm="rk",
                           ensemble_mth=c("mean",
                                          "top_k",
                                          "exp_k",
                                          "wt_mean",
                                          "wt_top_k",
                                          "wt_exp_k"),
                           wt=rep(1,length(rank_lst)),
                           k=NULL){
  
  #=== collect features and their rankings from different model
  p<-length(rank_lst) #number of different rankings
  varlst<-c()
  rk_stack<-c()
  for(i in 1:p){
    #subset and reorder columns so that var_colnm is the 1st column
    varlst_i<-rank_lst[[i]][,c(var_colnm,rank_colnm)]
    #rename column for easy referring
    varlst_i<-setNames(varlst_i,c("var","rk"))
    #collect all distinct features
    varlst<-unique(c(varlst,varlst_i$var))
    #stack feature rankings
    rk_stack %<>% 
      bind_rows(varlst_i %>% mutate(mod=paste0("m",i))) #add model index
  }
  n<-length(varlst) #number of distinct features ever got selected
  
  #=== impute rankings for un-selected features
  rank_df<-rk_stack %>%
    spread(mod,rk) %>%
    replace(is.na(.),Inf) %>%
    gather(mod,rk,-var) %>%
    group_by(mod) %>%
    mutate(rk=dense_rank(rk)) %>%
    # mutate(rk=rank(rk,ties.method = "min")) %>% # take the lowest ranking when ties occur
    # mutate(rk=rank(rk,ties.method = "max")) %>% # take the highest ranking when ties occur
    # mutate(rk=rank(rk,ties.method = "average")) %>% # take the average ranking when ties occur
    # mutate(rk=rank(rk,ties.method = "random")) %>% # randomly assign orders when ties occur
    ungroup 
  
  #==== rank transformation
  if(grepl("_k",ensemble_mth)){
    if(is.null(k)){
      stop("need to specify k, i.e. number of top features of interest!")
    }
    
    if(grepl("top_k",ensemble_mth)){
      rank_df %<>%
        mutate(rk=as.numeric(rk<=k)) %>%
        spread(var,rk)
    }else if(grepl("exp_k",ensemble_mth)){
      rank_df %<>%
        mutate(rk=(exp(-rk/k))) %>%
        spread(var,rk)
    }else{
      stop("current ensemble method is not supported!")
    }
  }
  
  if(grepl("wt$",ensemble_mth)){
    #weighted ensemble schemes
    if(p!=length(wt)){
      stop("length of wt should be the same as nrow(rank_mt)!")
    }
    
    if(all(wt==1)){
      warning("all weights in wt are 1!")
    }
    
    #broadcast the normalized weight vector
    rank_df<-data.frame(mod=rank_df[,1],
                        as.matrix(rank_df[,-1]) * wt/sum(wt))
  }
  
  #=== rank aggregation
  rank_df_pivot<-rank_df %>%
    gather(var,rk,-mod) %>%
    group_by(var) %>%
    summarise(agg_rk=mean(rk),
              .groups="drop") %>%
    ungroup %>%
    mutate(agg_rk=dense_rank(agg_rk))
  
  return(rank_df_pivot)
}

