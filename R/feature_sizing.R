###################################################################
# Feature Sizing Module                                           #
# Parameters:                                                     #
# -- X: a data.frame of predictors                                #
# -- y: a vector of real labels                                   #
# -- rank_lst:feature ranking list                                #
#    - list: each list can be a data.frame or matrix              #
#       with two columns (var_colnm, rank_colnm)                  #
# -- var_colnm: column name for variable in each list             #
# -- rank_colnm: column name for rank in each list                #
# -- model: type of base model ("xgb","glmnet","nnet","rf","svm") #
# -- verb: verbose, report progress                               #
###################################################################

source("./R/util.R")
source("./R/feature_ranking.R")

feature_sizing<-function(df,             #a data.frame: require "id","y" columns, while all the rest are predictors
                         var_rk,
                         var_colnm="var",
                         rank_colnm="rk",
                         model=list("xgb",
                                    "glmnet",
                                    "nnet",
                                    "rf",
                                    "svm"), #allow external models?
                         verb=T
                         ){
  start_fs<-Sys.time() #benchmark

  #initialize loop parameters
  a<-feat_sel_k_low
  d<-feat_sel_k_up
  b<-a+0.1
  c<-d-0.1
  auc_inc<--Inf
  inc_p<-1
  opt_size<-1
  global_min<-d
  ROC_opt_update<-0
  
  #initialize search path
  feat_num<-list()
  track_path<-c()
  
  #initialize AUROC as if it is by random chance
  ROC_obj_new<-pROC::roc(y,sample(c(0,1),length(y),replace=T),direction="<") 
  ROC_obj_opt<-ROC_obj_new
  
  #convert data to long-skinny format
  df_long<-df %>%
    gather(var,val,-id)
  
  #rename column for easy referring
  var_rk<-setNames(var_rk,c("var","rk"))
  
  #main body of golden-section-search algorithm
  while(a < (d-1) && b < (c-1)){
    start_g<-Sys.time() #benchmark
    if(verb) cat("start golden-section search \n")
    
    #track the approximation path
    track_path<-rbind(track_path,cbind(a=a,d=d))
    
    #update golden-section points: b,c
    b<-floor(d-(d-a)/phi)
    c<-floor(a+(d-a)/phi)
    bounds<-c(a,b,c,d)
    
    for(k in seq_along(bounds)){
      feat_sel_k<-bounds[k]
      
      if(!is.null(feat_num[[paste0("size_",feat_sel_k)]])){
        if(verb) cat("...the case of keeping",feat_sel_k,"features has already been saved. \n")
      }else{
        start_k<-Sys.time()
        if(verb) cat("...keep",feat_sel_k,"features and retrain the model. \n")
        
        #subset features
        df_subset<-df_long %>% 
          semi_join(var_rk %>%
                      arrange(rk) %>% 
                      dplyr::select(var,rk) %>% 
                      dplyr::slice(1:feat_sel_k), 
                    by=c("var")) %>% 
          spread(var,val)
        
        feat_rk_out<-feature_ranking(data=df_subset,
                                     ...)
        
        #standard output of feature_ranking():
        # - hyper_param
        # - model
        # - pred_df
        # - var_rk
        # - benchmark
        
        # evaluate auc improvement
        ROC_obj_new<-pROC::auc(feat_rk_out$pred_df$real,
                               feat_rk_out$pred_df$pred)

        #need to update opt?
        if(ROC_obj_new$auc > ROC_obj_opt$auc){
          ROC_obj_opt<-ROC_obj_new
          opt_size<-feat_sel_k
          ROC_opt_update<-ROC_opt_update+1
        }
        
        # memorize this scenario, in case a callback
        feat_num[[paste0("size_",feat_sel_k)]]<-list(roc_obj=ROC_obj_new,
                                                     ROC_opt_update=ROC_opt_update)
        
        #-----------------------------------------------benchmark------------------------------------------------------------#
        time_perf_i_nm<-c(time_perf_i_nm,paste0("finish_eval_feature_num@",feat_sel_k))
        time_perf_i<-c(time_perf_i,paste0(Sys.time()-start_k,units(Sys.time()-start_k)))
        if(verb) cat("...finish modeling on",feat_sel_k,"features for",fs_mth[fs],"in",time_perf_i[length(time_perf_i)],"\n")
        #-----------------------------------------------benchmark------------------------------------------------------------#
      }
    } #end for loop
      
    #compare b(2),c(3)
    auc_b_c<-feat_num[[paste0("size_",bounds[2])]]$roc_obj$auc-feat_num[[paste0("size_",bounds[3])]]$roc_obj$auc
    aucp_b_c<-pROC::roc.test(feat_num[[paste0("size_",bounds[2])]]$roc_obj,feat_num[[paste0("size_",bounds[3])]]$roc_obj,method='delong')$p.value
    
    #compare b(2) with opt
    auc_b_opt<-feat_num[[paste0("size_",bounds[2])]]$roc_obj$auc-ROC_obj_opt$auc
    aucp_b_opt<-pROC::roc.test(feat_num[[paste0("size_",bounds[2])]]$roc_obj,ROC_obj_opt,method='delong')$p.value
    
    #compare c with opt
    auc_c_opt<-feat_num[[paste0("size_",bounds[3])]]$roc_obj$auc-ROC_obj_opt$auc
    aucp_c_opt<-pROC::roc.test(feat_num[[paste0("size_",bounds[3])]]$roc_obj,ROC_obj_opt,method='delong')$p.value
    
    #update a,b,c,d
    if((max(aucp_b_opt,aucp_c_opt) <= inc_tol_p &&
        max(auc_b_opt,auc_c_opt) < 0)){
      a<-b
      d<-opt_size
      local_min<<-opt_size
    }else if(aucp_b_c > inc_tol_p ||
             auc_b_c > 0){
      d<-c
      local_min<<-b
    }else if(aucp_c_opt > inc_tol_p ||
             auc_c_opt >= 0){
      a<-b
      local_min<<-c
    }else{
      stop("conditions are not exhaustive!")
    }
    
    #update global_min?
    if(local_min < global_min){
      global_min <- local_min
      min_model <<- feat_rk_out$model
    }
    
    #--------------------------------------benchmark---------------------------------------------------------#
    time_perf_i_nm<-c(time_perf_i_nm,paste0("shrink_search_interval_to_",a,"_",d))
    time_perf_i<-c(time_perf_i,paste0(Sys.time()-start_g,units(Sys.time()-start_g)))
    if(verb) cat(fs_mth[fs],":shrink interval to[",a,",",d,"] \n")
    #--------------------------------------benchmark---------------------------------------------------------#
  
  } #end while loop
  feat_num$opt_model<-feat_rk_out #only record the model with optimal feature size
  feat_num$track_path<-track_path #record the search track
  feat_num$benchmark<-data.frame(benchmark_task=time_perf_i_nm,
                                 benchmark_time=time_perf_i,
                                 stringsAsFactors = F)
  return(feat_num)
}
