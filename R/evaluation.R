##########################################################
# Evaluation Module                                      #
# - eval_stability(): evaluate feature ranking stability #
#   parameters:                                          #
#   -- rank_lst:feature ranking list                     #
#    - list: each list can be a data.frame or matrix     #
#       with two columns (var_colnm, rank_colnm)         #
#   -- var_colnm: column name for variable in each list  #
#   -- rank_colnm: column name for rank in each list     #
#                                                        #
# - eval_accuracy(): evaluate prediction accuracy        #
#   parameters:                                          #
#   -- real: real label from testing set                 #
#   -- pred: predicted label from testing set            #
#   -- keep_all_cutoffs: the full ROC curve is preserved #
#                        if set to "T"                   #  
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


eval_accuracy<-function(pred,
                        real,
                        keep_all_cutoffs=F){
  # various performace table
  pred_obj<-ROCR::prediction(pred,real)
  
  prc<-performance(pred_obj,"prec","rec")
  roc<-performance(pred_obj,"sens","spec")
  nppv<-performance(pred_obj,"ppv","npv")
  pcfall<-performance(pred_obj,"pcfall")
  acc<-performance(pred_obj,"acc")
  fscore<-performance(pred_obj,"f")
  mcc<-performance(pred_obj,"phi")
  
  perf_at<-data.frame(cutoff=prc@alpha.values[[1]],
                      prec=prc@y.values[[1]],
                      rec_sens=prc@x.values[[1]],
                      stringsAsFactors = F) %>% 
    arrange(cutoff) %>%
    left_join(data.frame(cutoff=nppv@alpha.values[[1]],
                         ppv=nppv@y.values[[1]],
                         npv=nppv@x.values[[1]],
                         stringsAsFactors = F),
              by="cutoff") %>%
    dplyr::mutate(prec_rec_dist=abs(prec-rec_sens)) %>%
    left_join(data.frame(cutoff=fscore@x.values[[1]],
                         fscore=fscore@y.values[[1]],
                         stringsAsFactors = F),
              by="cutoff") %>%
    left_join(data.frame(cutoff=roc@alpha.values[[1]],
                         spec=roc@x.values[[1]],
                         stringsAsFactors = F),
              by="cutoff") %>%
    dplyr::mutate(Euclid_meas=sqrt((1-rec_sens)^2+(0-(1-spec))^2),
                  Youden_meas=rec_sens+spec-1) %>%
    left_join(data.frame(cutoff=pcfall@x.values[[1]],
                         pcfall=pcfall@y.values[[1]],
                         stringsAsFactors = F),
              by="cutoff") %>%
    left_join(data.frame(cutoff=acc@x.values[[1]],
                         acc=acc@y.values[[1]],
                         stringsAsFactors = F),
              by="cutoff") %>%
    left_join(data.frame(cutoff=mcc@x.values[[1]],
                         mcc=mcc@y.values[[1]],
                         stringsAsFactors = F),
              by="cutoff") %>%
    filter(prec > 0 & rec_sens > 0 & spec > 0) %>%
    group_by(cutoff) %>%
    dplyr::mutate(size=n()) %>%
    ungroup
  
  # performance summary
  lab1<-pred[real==1]
  lab0<-pred[real==0]
  pr<-pr.curve(scores.class0 = lab1,
               scores.class1 = lab0,curve=F)
  roc_ci<-pROC::ci.auc(real,pred)
  
  perf_summ<-data.frame(overall_meas=c("roauc_low",
                                       "roauc",
                                       "roauc_up",
                                       "opt_thresh",
                                       "opt_sens",
                                       "opt_spec",
                                       "opt_ppv",
                                       "opt_npv",
                                       "prauc1",
                                       "prauc2",
                                       "opt_prec",
                                       "opt_rec",
                                       "opt_fscore"),
                        meas_val=c(roc_ci[[1]],
                                   roc_ci[[2]],
                                   roc_ci[[3]],
                                   perf_at$cutoff[which.min(perf_at$Euclid_meas)],
                                   perf_at$rec_sens[which.min(perf_at$Euclid_meas)],
                                   perf_at$spec[which.min(perf_at$Euclid_meas)],
                                   perf_at$ppv[which.min(perf_at$Euclid_meas)],
                                   perf_at$npv[which.min(perf_at$Euclid_meas)],
                                   pr$auc.integral,
                                   pr$auc.davis.goadrich,
                                   perf_at$prec[which.min(perf_at$prec_rec_dist)],
                                   perf_at$rec_sens[which.min(perf_at$prec_rec_dist)],
                                   perf_at$fscore[which.min(perf_at$prec_rec_dist)]),
                        stringsAsFactors = F) %>%
    bind_rows(perf_at %>% 
                dplyr::summarize(prec_m=mean(prec,na.rm=T),
                                 sens_m=mean(rec_sens,na.rm=T),
                                 spec_m=mean(spec,na.rm=T),
                                 ppv_m=mean(ppv,na.rm=T),
                                 npv_m=mean(npv,na.rm=T),
                                 acc_m=mean(acc,na.rm=T),
                                 fscore_m=mean(fscore,na.rm=T),
                                 mcc_m=mean(mcc,na.rm=T),
                                 .groups="drop") %>%
                gather(overall_meas,meas_val))
  
  out<-list(perf_summ=perf_summ)
  if(keep_all_cutoffs){
    out$perf_at<-perf_at
  }
  
  return(out)
}


