##########################################################
# Evaluation Module                                      #
# - eval_stability(): evaluate feature ranking stability #
# - parameters:                                          #
# -- rank_lst:feature ranking list                       #
#    - list: each list can be a data.frame or matrix     #
#       with two columns (var_colnm, rank_colnm)         #
# -- var_colnm: column name for variable in each list    #
# -- rank_colnm: column name for rank in each list       #
#                                                        #
# - eval_accuracy(): evaluate prediction accuracy        #
##########################################################

eval_stability<-function(rank_lst,
                         var_colnm="var",
                         rank_colnm="rk",
                         metric_type=c("kuncheva",
                                       "wci",
                                       "wci_rel"),
                         f=1000,
                         d=200){
  
  K<-length(rank_lst)
  rk_stack<-c()
  for(i in 1:K){
    #subset and reorder columns so that var_colnm is the 1st column
    varlst_i<-rank_lst[[i]][,c(var_colnm,rank_colnm)]
    #rename column for easy referring
    varlst_i<-setNames(varlst_i,c("var","rk"))
    #collect all distinct features
    varlst<-unique(c(varlst,varlst_i$var))
    #stack feature rankings
    rk_stack %<>% 
      bind_rows(varlst_i %>% 
                  mutate(mod_idx=i) %>%
                  mutate(p=nrow(varlst_i))) #add model index
  }
  
  if(metric_type="kuncheva"){
    if(is.null(f)|is.null(d)){
      stop("need to specify the total feature size, f; and a fixed subset size, d")
    }
    
    #kalousis and kuncheva index
    pairwise_fset<-c()
    for(i in 1:(K-1)){
      for(j in (i+1):K){
        ki_new<-rk_stack %>%
          filter(mod_idx %in% c(i,j)) %>%
          filter(rk<=pmin(p,d)) %>%
          dplyr::select(var,mod) %>%
          mutate(pair_idx=paste0(i,"_",j)) %>%
          mutate(mod_idx=ifelse(mod_idx==i,"ki","kj"),ref=1) %>%
          unique %>% spread(mod_idx,ref,fill=0) %>%
          mutate(inter=ki*kj,union=((ki+kj)>=1)*1)
        
        pairwise_fset %<>%
          bind_rows(ki_new)
      }
    }
    
    stb_idx<-pairwise_fset %>%
      group_by(pair_idx) %>%
      summarize(inter_cnt=sum(inter),
                union_cnt=sum(union)) %>%
      ungroup %>%
      summarize(ki_unscale=sum((inter_cnt*f-d^2)/(d*(f-d)))*(2/(K*(K-1)))) %>%
      mutate(ki=(ki_unscale-(-1))/(1-(-1))) #rescale
  }
  
  else if(metric_type="wci"){
    #weighted consistency index
    stb_idx<-rk_stack %>%
      group_by(var) %>% 
      dplyr::summarize(phi_f=n()) %>%
      ungroup %>% 
      dplyr::mutate(N=sum(phi_f)) %>%
      dplyr::mutate(phi_f_wt=(phi_f/N)*((phi_f-1)/(K-1))) %>%
      group_by(N) %>%
      dplyr::summarize(cw=sum(phi_f_wt),
                       C=length(unique(var))) %>%
      ungroup %>%
      mutate(D=N %% C,
             H=N %% K) %>%
      mutate(cw_min=(N^2-C*(N-D)-D^2)/(C*N*(K-1)),
             cw_max=(H^2+N*(K-1)-H*K)/(N*(K-1))) %>%
      mutate(cw_rel=(cw-cw_min)/(cw_max-cw_min)) %>%
      select(cw,cw_rel) 
  }
  
  else if(metric_type="wci_rel"){
    #(relative) weighted consistency index
    stb_idx<-rk_stack %>%
      group_by(p,var) %>% 
      dplyr::summarize(phi_f=n()) %>%
      ungroup %>% 
      group_by(p) %>%
      dplyr::mutate(N=sum(phi_f)) %>%
      ungroup %>%
      dplyr::mutate(phi_f_wt=(phi_f/N)*((phi_f-1)/(K-1))) %>%
      group_by(p,N) %>%
      dplyr::summarize(cw=sum(phi_f_wt),
                       C=length(unique(Feature))) %>%
      ungroup %>%
      mutate(D=N %% C,
             H=N %% K) %>%
      mutate(cw_min=(N^2-C*(N-D)-D^2)/(C*N*(K-1)),
             cw_max=(H^2+N*(K-1)-H*K)/(N*(K-1))) %>%
      mutate(cw_rel=(cw-cw_min)/(cw_max-cw_min)) %>%
      select(p,cw,cw_rel) 
  }
  
  else{
    stop("stability index calculation method is not supported!")
  }
  
  return(stb_idx)
}


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

