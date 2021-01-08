##########################################################
# Evaluation Module                                      #
# Parameters:                                            #
# -- rank_lst:feature ranking list                       #
#    - list: each list can be a data.frame or matrix     #
#       with two columns (var_colnm, rank_colnm)         #
# -- var_colnm: column name for variable in each list    #
# -- rank_colnm: column name for rank in each list       #
##########################################################

eval_accuracy<-function(){
  ####################################################################
  #calculate accuracy for each step on the search path
  start<-Sys.time()
  accuracy_stack<-c()
  opt_fs<-c()
  for(i in 1:10){
    start_i<-Sys.time()
    
    for(fs in seq_along(fs_mth)){
      fs_i<-fs_summary[[paste0("resample",i,"@",fs_mth[fs])]]
      for(j in seq_len(length(fs_i)-2)){
        new_stack<-fs_i[[j]][["model_summary"]] %>%
          dplyr::select(PATIENT_NUM,part73,real,pred,fs_num,fs) %>%
          filter(part73!="T") %>% #remove training points
          group_by(fs,fs_num,part73) %>%
          dplyr::summarize(auc=auc(real,pred),
                           thre_ss=paste(coords(roc(real,pred),x="best")[1:3],collapse=",")) %>%
          separate(thre_ss,c("thresh","spec","sens"),sep=",",remove=F) %>% 
          select(-thre_ss) %>% mutate(auc=as.numeric(auc),
                                      thresh=as.numeric(thresh),
                                      sens=as.numeric(sens),
                                      spec=as.numeric(spec)) %>% 
          ungroup %>% mutate(resample = i)
        
        accuracy_stack %<>%
          bind_rows(new_stack)
      }
      
      opt_fs %<>%
        bind_rows(new_stack %>% 
                    dplyr::select(fs_num,fs,resample) %>% unique)
    }
    lapse_i<-Sys.time()-start_i
    cat("finish evaluating accuracy for resample",i,"in",lapse_i,units(lapse_i),"\n")
  }
  Sys.time()-start 
}

eval_stability<-function(){
  ##############################################################################################################
  #calculate stability for a system (adaptive feature size)
  K<-resamples
  fs_mth_elig<-unique(accuracy_stack$fs)
  
  feat_stack<-c()
  for(i in 1:5){
    opt_fs_i<-opt_fs %>% 
      filter(resample==i) %>% select(fs,fs_num)
    
    feat_stack %<>%
      bind_rows(feature_rk[[paste0("resample",i)]] %>%
                  gather(fs_mth,adj_rank,-Feature) %>%
                  dplyr::filter(!is.na(adj_rank) & (fs_mth %in% fs_mth_elig)) %>%
                  left_join(opt_fs_i,by=c("fs_mth"="fs")) %>%
                  group_by(fs_mth) %>%
                  dplyr::mutate(rank=rank(adj_rank,ties.method="random")) %>%
                  filter(rank <= fs_num) %>% 
                  dplyr::select(-rank,-adj_rank) %>%
                  dplyr::mutate(resample=i))
  }
  
  #######################################
  #(relative) weighted consistency index#
  #######################################
  cw_idx<-feat_stack %>%
    group_by(fs_mth,Feature) %>% 
    dplyr::summarize(phi_f=n()) %>%
    ungroup %>% 
    group_by(fs_mth) %>%
    dplyr::mutate(N=sum(phi_f)) %>%
    ungroup %>%
    dplyr::mutate(phi_f_wt=(phi_f/N)*((phi_f-1)/(K-1))) %>%
    group_by(fs_mth,N) %>%
    dplyr::summarize(cw=sum(phi_f_wt),
                     C=length(unique(Feature))) %>%
    ungroup %>%
    mutate(D=N %% C,
           H=N %% K) %>%
    mutate(cw_min=(N^2-C*(N-D)-D^2)/(C*N*(K-1)),
           cw_max=(H^2+N*(K-1)-H*K)/(N*(K-1))) %>%
    mutate(cw_rel=(cw-cw_min)/(cw_max-cw_min)) %>%
    select(fs_mth,cw,cw_rel) %>%
    ungroup %>%
    mutate(fs_mth = recode(fs_mth,
                           single_rank="single_run_rank",
                           mean_rank="ensemble1_mean_rank",
                           wt_mean_rank="ensemble2_mean_rank_wt",
                           # top_50="ensemble2_top_50",
                           # wt_top_50="ensemble2_top_50_wt",
                           top_100="ensemble3_top_100",
                           wt_top_100="ensemble4_top_100_wt",
                           # top_50_exp="ensemble4_top_50_exp",
                           # wt_top_50_exp="ensemble4_top_50_exp_wt",
                           top_100_exp="ensemble5_top_100_exp",
                           wt_top_100_exp="ensemble6_top_100_exp_wt")) #recode fs_mth for better legend in plots
  
  #save tables
  stable_stack_adapt<-list(cw_idx=cw_idx)
  ######################################################################################################################################
  
  ######################################################################################################################################
  #stability vs accuracy at optimal feature size
  accu_stab_opt<-stable_stack_adapt$cw_idx %>%
    left_join(accuracy_summ$accuracy_optsize %>%
                filter(part73=="V") %>%
                select(resample,fs,fs_num,auc) %>% 
                group_by(fs,fs_num) %>%
                filter(auc==max(auc)) %>%
                unique,
              by=c("fs_mth"="fs")) %>%
    mutate(euclid=sqrt(cw_rel^2+(auc-0.5)^2),
           geo_root=cw_rel*(auc-0.5))
  
  save(accu_stab_opt,file=paste0(model_type,"accuracy_stability_opt.Rdata"))
  load(paste0(model_type,"accuracy_stability_opt.Rdata"))
  ##############################################################################################################
  
  ##############################################################################################################
  #calculate stability for a system (assuming fixed feature size)
  K<-resamples
  fs_mth_elig<-unique(accuracy_stack$fs)
  brks<-c(20,50,100,150,200,250,300,350,400,450,500)
  
  feat_stack<-c()
  for(i in 1:resamples){
    for(j in brks) {
      feat_stack %<>%
        bind_rows(feature_rk[[paste0("resample",i)]] %>%
                    gather(fs_mth,adj_rank,-Feature) %>%
                    filter(!is.na(adj_rank) & (fs_mth %in% fs_mth_elig)) %>%
                    group_by(fs_mth) %>%
                    mutate(rank=rank(adj_rank,ties.method="random")) %>%
                    mutate(fs_num=cut(rank,breaks=c(1,j),labels=j,include.lowest=T)) %>%
                    filter(!is.na(fs_num)) %>% dplyr::select(-rank,-adj_rank) %>% unique %>%
                    mutate(fs_num=as.numeric(as.character(fs_num)),resample=i))
    }
    
  }
  
  #######################################
  #(relative) weighted consistency index#
  #######################################
  cw_idx<-feat_stack %>%
    group_by(fs_mth,fs_num,Feature) %>% 
    dplyr::summarize(phi_f=n()) %>%
    ungroup %>% 
    group_by(fs_mth,fs_num) %>%
    dplyr::mutate(N=sum(phi_f)) %>%
    ungroup %>%
    dplyr::mutate(phi_f_wt=(phi_f/N)*((phi_f-1)/(K-1))) %>%
    group_by(fs_mth,fs_num,N) %>%
    dplyr::summarize(cw=sum(phi_f_wt),
                     C=length(unique(Feature))) %>%
    ungroup %>%
    mutate(D=N %% C,
           H=N %% K) %>%
    mutate(cw_min=(N^2-C*(N-D)-D^2)/(C*N*(K-1)),
           cw_max=(H^2+N*(K-1)-H*K)/(N*(K-1))) %>%
    mutate(cw_rel=(cw-cw_min)/(cw_max-cw_min)) %>%
    select(fs_mth,fs_num,cw,cw_rel) %>%
    ungroup %>%
    mutate(fs_mth = recode(fs_mth,
                           single_rank="single_run_rank",
                           mean_rank="ensemble1_mean_rank",
                           wt_mean_rank="ensemble2_mean_rank_wt",
                           # top_50="ensemble2_top_50",
                           # wt_top_50="ensemble2_top_50_wt",
                           top_100="ensemble3_top_100",
                           wt_top_100="ensemble4_top_100_wt",
                           # top_50_exp="ensemble4_top_50_exp",
                           # wt_top_50_exp="ensemble4_top_50_exp_wt",
                           top_100_exp="ensemble5_top_100_exp",
                           wt_top_100_exp="ensemble6_top_100_exp_wt")) #recode fs_mth for better legend in plots
  
  
  #########################################
  #generalized Kalousis and Kuncheva index#
  #########################################
  tot_feat_cnt<-600*10
  
  pairwise_fset<-c()
  for(i in 1:(K-1)){
    for(j in (i+1):K){
      ki_new<-feat_stack %>%
        filter(resample %in% c(i,j)) %>%
        dplyr::select(Feature,fs_mth,fs_num,resample) %>%
        mutate(pair_idx=paste0(i,"_",j)) %>%
        mutate(resample=ifelse(resample==i,"ki","kj"),ref=1) %>%
        unique %>% spread(resample,ref,fill=0) %>%
        mutate(inter=ki*kj,union=((ki+kj)>=1)*1)
      
      pairwise_fset %<>%
        bind_rows(ki_new)
    }
  }
  
  kch_idx<-pairwise_fset %>%
    group_by(fs_mth,fs_num,pair_idx) %>%
    summarize(inter_cnt=sum(inter),
              union_cnt=sum(union)) %>%
    group_by(fs_mth,fs_num) %>%
    summarize(ki_unscale=sum((inter_cnt*tot_feat_cnt-fs_num^2)/(fs_num*(tot_feat_cnt-fs_num)))*(2/(K*(K-1)))) %>%
    mutate(ki=(ki_unscale-(-1))/(1-(-1))) %>% #scale index to 0,1 range
    ungroup %>%
    mutate(fs_mth = recode(fs_mth,
                           single_rank="single_run_rank",
                           mean_rank="ensemble1_mean_rank",
                           wt_mean_rank="ensemble2_mean_rank_wt",
                           # top_50="ensemble2_top_50",
                           # wt_top_50="ensemble2_top_50_wt",
                           top_100="ensemble3_top_100",
                           wt_top_100="ensemble4_top_100_wt",
                           # top_50_exp="ensemble4_top_50_exp",
                           # wt_top_50_exp="ensemble4_top_50_exp_wt",
                           top_100_exp="ensemble5_top_100_exp",
                           wt_top_100_exp="ensemble6_top_100_exp_wt")) #recode fs_mth for better legend in plots
  
}

