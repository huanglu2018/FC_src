options(stringsAsFactors = F)
library(plyr)
library(tidyverse)
library(doParallel)
library(magrittr)
library(data.table)
library(parallel)
library(ggplot2)
library(Hmisc)



# Yuan et.al. signature 6 genes

get_score_M1=function(df){
  M_gene=c("AHCYL2","LAMP2","SPRY1","SERPINA7","FGGY","YBX1P4")
  # id_transfer(M_gene,"ENSEMBL","SYMBOL")
  tumor_gene_df=df[,grep("^nat__",colnames(df),value = T,invert = T)] %>% 
    set_colnames(gsub("tumor__","",colnames(.)))
  MISSGENE=setdiff(M_gene,colnames(tumor_gene_df))
  if(length(MISSGENE)==0){
    model_gene_df=tumor_gene_df[,c("TIME_TO_EVENT","EVENT",M_gene)]
    model_gene_df %>% 
      rownames_to_column(var="ind") %>% 
      mutate(risk=0.51*AHCYL2+0.54*LAMP2+0.36*SPRY1+0.33*SERPINA7+0.33*FGGY+0.18*YBX1P4+0.001) %>% 
      dplyr::select(c("TIME_TO_EVENT","EVENT","risk","ind")) %>% 
      column_to_rownames(var="ind")
  }else{
    cat(paste0("genes lost: ",paste0(MISSGENE,collapse = ",")))
  }
  
}
M1_rs=get_score_M1(vali_es_df)



# r Kim et.al. signature 65 genes

kim_model=read.csv(kim_model_file)[,1:2] %>% set_colnames(c("gene","coef"))

get_score_M2=function(df,kim_model,normalize=F){
  tumor_gene_df=df[,grep("^nat__",colnames(df),value = T,invert = T)] %>% 
    set_colnames(gsub("tumor__","",colnames(.)))
  MISSGENE=setdiff(kim_model$Gene,colnames(tumor_gene_df))
  
  if(length(MISSGENE)==0){
    model_gene_df=tumor_gene_df[,c("TIME_TO_EVENT","EVENT",kim_model$gene)] %>% na.omit()
    model_gene_df[,"risk"]=NA 
    for (i in 1:nrow(model_gene_df)) { 
      sub_risk_info=model_gene_df[i,kim_model$gene] %>% 
        t() %>% 
		as.data.frame() %>% 
		set_colnames("ind") %>% 
        merge(.,kim_model,by.x=0,by.y="gene") %>% 
		mutate(sub_risk=ind*coef) 
      if(normalize==F){ 
        model_gene_df[i,"risk"]=sum(sub_risk_info$sub_risk) 
      }else{ 
        model_gene_df[i,"risk"]=sum(sub_risk_info$sub_risk)/sum(sub_risk_info$ind) 
      }
      # print(sum(sub_risk_info$sub_risk))
    }
    model_gene_df %>% dplyr::select(c("TIME_TO_EVENT","EVENT","risk"))
    
  }else{
    cat(paste0("genes lost: ",paste0(MISSGENE,collapse = ",")))
  }
}
M2_rs=get_score_M2(vali_es_df,kim_model)


# jiang et.al. 5 gene signature

jiang_model=read.csv(jiang_model_file)[,1:2] %>% set_colnames(c("coef","gene"))

get_score_M3=function(df,kim_model,normalize=F){
  tumor_gene_df=df[,grep("^nat__",colnames(df),value = T,invert = T)] %>% 
    set_colnames(gsub("tumor__","",colnames(.)))
  MISSGENE=setdiff(kim_model$Gene,colnames(tumor_gene_df))
  
  if(length(MISSGENE)==0){
    model_gene_df=tumor_gene_df[,c("TIME_TO_EVENT","EVENT",kim_model$gene)] %>% na.omit()
    # model_gene_df %>% mutate(risk=0.51*AHCYL2+0.54*LAMP2+0.36*SPRY1+0.33*SERPINA7+0.33*FGGY+0.18*YBX1P4+0.001) %>% 
    #   select(c("TIME_TO_EVENT","EVENT","risk"))
    model_gene_df[,"risk"]=NA
    for (i in 1:nrow(model_gene_df)) {
      sub_risk_info=model_gene_df[i,kim_model$gene] %>% 
	    t() %>% 
		as.data.frame() %>% 
		set_colnames("ind") %>% 
        merge(.,kim_model,by.x=0,by.y="gene") %>% 
		mutate(sub_risk=ind*coef)
      if(normalize==F){
        model_gene_df[i,"risk"]=sum(sub_risk_info$sub_risk)
      }else{
        model_gene_df[i,"risk"]=sum(sub_risk_info$sub_risk)/sum(sub_risk_info$ind)
      }
      # print(sum(sub_risk_info$sub_risk))
    }
    model_gene_df %>% dplyr::select(c("TIME_TO_EVENT","EVENT","risk"))
    
  }else{
    cat(paste0("genes lost: ",paste0(MISSGENE,collapse = ",")))
  }
}
M3_rs=get_score_M3(vali_es_df,jiang_model)


# kong 3 gene signature

kong_model=read.table(kong_model_file,header = T)[,1:2] %>%
  set_colnames(c("coef","gene")) %>% mutate(coef=as.numeric(coef))

get_score_M4=function(df,kong_model,normalize=F){
  tumor_gene_df=df[,grep("^nat__",colnames(df),value = T,invert = T)] %>% 
    set_colnames(gsub("tumor__","",colnames(.)))
  MISSGENE=setdiff(kong_model$Gene,colnames(tumor_gene_df))
  
  if(length(MISSGENE)==0){
    model_gene_df=tumor_gene_df[,c("TIME_TO_EVENT","EVENT",kong_model$gene)] %>% na.omit()
    # model_gene_df %>% mutate(risk=0.51*AHCYL2+0.54*LAMP2+0.36*SPRY1+0.33*SERPINA7+0.33*FGGY+0.18*YBX1P4+0.001) %>% 
    #   select(c("TIME_TO_EVENT","EVENT","risk"))
    model_gene_df[,"risk"]=NA
    for (i in 1:nrow(model_gene_df)) {
      sub_risk_info=model_gene_df[i,kong_model$gene] %>% 
	    t() %>% 
		as.data.frame() %>% 
		set_colnames("ind") %>% 
        merge(.,kong_model,by.x=0,by.y="gene") %>% 
		mutate(sub_risk=ind*coef)
      if(normalize==F){
        model_gene_df[i,"risk"]=sum(sub_risk_info$sub_risk)
      }else{
        model_gene_df[i,"risk"]=sum(sub_risk_info$sub_risk)/sum(sub_risk_info$ind)
      }
      # print(sum(sub_risk_info$sub_risk))
    }
    model_gene_df %>% dplyr::select(c("TIME_TO_EVENT","EVENT","risk"))
    
  }else{
    cat(paste0("genes lost: ",paste0(MISSGENE,collapse = ",")))
  }
}
M4_rs=get_score_M4(vali_es_df,kong_model)



# liu 7 gene signature

liu_model=read.table(liu_model_file,header = T)[,1:2] %>% 
  set_colnames(c("gene","coef")) %>%
  mutate(coef=as.numeric(coef))

get_score_M5=function(df,liu_model,normalize=F){
  tumor_gene_df=df[,grep("^nat__",colnames(df),value = T,invert = T)] %>% 
    set_colnames(gsub("tumor__","",colnames(.)))
  MISSGENE=setdiff(liu_model$Gene,colnames(tumor_gene_df))
  
  if(length(MISSGENE)==0){
    model_gene_df=tumor_gene_df[,c("TIME_TO_EVENT","EVENT",liu_model$gene)] %>% na.omit()
    model_gene_df[,"risk"]=NA
    for (i in 1:nrow(model_gene_df)) {
      sub_risk_info=model_gene_df[i,liu_model$gene] %>% 
	    t() %>%
	    as.data.frame() %>% 
	    set_colnames("ind") %>% 
        merge(.,liu_model,by.x=0,by.y="gene") %>% 
		mutate(sub_risk=ind*coef)
      if(normalize==F){
        model_gene_df[i,"risk"]=sum(sub_risk_info$sub_risk)
      }else{
        model_gene_df[i,"risk"]=sum(sub_risk_info$sub_risk)/sum(sub_risk_info$ind)
      }
      # print(sum(sub_risk_info$sub_risk))
    }
    model_gene_df %>% dplyr::select(c("TIME_TO_EVENT","EVENT","risk"))
    
  }else{
    cat(paste0("genes lost: ",paste0(MISSGENE,collapse = ",")))
  }
}
M5_rs=get_score_M5(vali_es_df,liu_model)

model7_list=list(my=my_surv_df,
					M1=M1_rs,
					M2=M2_rs,
					M3=M3_rs,
					M4=M4_rs,
					M5=M5_rs,
					TG1=TG1_es,
					TG2=TG2_es,
					TG3=TG3_es,
					TG4=TG4_es)


my=model7_list$my %>% plyr::rename(c("GES_score"="my_risk"))
m1=model7_list$M1 %>% plyr::rename(c("risk"="m1_risk")) %>% .[,3,drop=F]
m2=model7_list$M2 %>% plyr::rename(c("risk"="m2_risk")) %>% .[,3,drop=F]
m3=model7_list$M3 %>% plyr::rename(c("risk"="m3_risk")) %>% .[,3,drop=F]
m4=model7_list$M4 %>% plyr::rename(c("risk"="m4_risk")) %>% .[,3,drop=F]
m5=model7_list$M5 %>% plyr::rename(c("risk"="m5_risk")) %>% .[,3,drop=F]
TG1=model7_list$TG1 %>% plyr::rename(c("risk"="TG1_risk")) %>% .[,3,drop=F]
TG2=model7_list$TG2 %>% plyr::rename(c("risk"="TG2_risk")) %>% .[,3,drop=F]
TG3=model7_list$TG3 %>% plyr::rename(c("risk"="TG3_risk")) %>% .[,3,drop=F]
TG4=model7_list$TG4 %>% plyr::rename(c("risk"="TG4_risk")) %>% .[,3,drop=F]

merged_risk=merge(my,m1,by=0) %>% 
  merge(.,m2,by.x=1,by.y=0) %>% 
  merge(.,m3,by.x=1,by.y=0) %>% 
  merge(.,m4,by.x=1,by.y=0) %>% 
  merge(.,m5,by.x=1,by.y=0) %>% 
  merge(.,TG1,by.x=1,by.y=0) %>% 
  merge(.,TG2,by.x=1,by.y=0) %>% 
  merge(.,TG3,by.x=1,by.y=0) %>% 
  merge(.,TG4,by.x=1,by.y=0)

library(prodlim)
library(survival)
# dat <- SimSurv(100)

# fit some candidate Cox models and compute the Kaplan-Meier estimate
Models <- list("my_model"=coxph(Surv(TIME_TO_EVENT,EVENT)~my_risk,data=merged_risk,x=TRUE,y=TRUE),
               "M1"=coxph(Surv(TIME_TO_EVENT,EVENT)~m1_risk,data=merged_risk,x=TRUE,y=TRUE),
               "M2"=coxph(Surv(TIME_TO_EVENT,EVENT)~m2_risk,data=merged_risk,x=TRUE,y=TRUE),
               "M3"=coxph(Surv(TIME_TO_EVENT,EVENT)~m3_risk,data=merged_risk,x=TRUE,y=TRUE),
               "M4"=coxph(Surv(TIME_TO_EVENT,EVENT)~m4_risk,data=merged_risk,x=TRUE,y=TRUE),
               "M5"=coxph(Surv(TIME_TO_EVENT,EVENT)~m5_risk,data=merged_risk,x=TRUE,y=TRUE),
               "TG1"=coxph(Surv(TIME_TO_EVENT,EVENT)~TG1_risk,data=merged_risk,x=TRUE,y=TRUE),
               "TG2"=coxph(Surv(TIME_TO_EVENT,EVENT)~TG2_risk,data=merged_risk,x=TRUE,y=TRUE),
               "TG3"=coxph(Surv(TIME_TO_EVENT,EVENT)~TG3_risk,data=merged_risk,x=TRUE,y=TRUE),
               "TG4"=coxph(Surv(TIME_TO_EVENT,EVENT)~TG4_risk,data=merged_risk,x=TRUE,y=TRUE)
)
# compute the apparent prediction error
if(CV==F){
PredError <- pec(object=Models,
                 formula=Surv(TIME_TO_EVENT,EVENT)~1,
                 data=merged_risk,
                 exact=F,
                 cens.model="marginal",
                 splitcMethod="none",
                 testTimes=c(12,24,36),
                 B=0,
                 verbose=TRUE)
}else{
  PredError <- pec(object=Models,
                   formula=Surv(TIME_TO_EVENT,EVENT)~1,
                   data=merged_risk,
                   exact=F,
                   cens.model="marginal",
                   splitMethod="Boot632plus",
                   B=100,
                   testTimes=c(12,24,36),
                   verbose=TRUE)
}

print(PredError,seq(12,36,12))
summary(PredError,times=seq(12,36,12))
plot(PredError,xlim=c(0,36),alpha=0.5)
}

