options(stringsAsFactors = F)
library(plyr)
library(tidyverse)
library(doParallel)
library(magrittr)
library(data.table)
library(parallel)
library(ggplot2)
library(Hmisc)
library(survivalROC)
library(timrROC)

each_yr_roc=function(roc_data,ifsmooth=F){
  if(ifsmooth==T){
    roc_data2=roc_data %>% set_colnames(c("TIME_TO_EVENT","EVENT","risk"))
    AUC_12=draw_ROC(12,roc_data2,colName[1]); par(new=TRUE)
    AUC_24=draw_ROC(24,roc_data2,colName[2]); par(new=TRUE)
    AUC_36=draw_ROC(36,roc_data2,colName[3]); par(new=TRUE)
    legend("bottomright", 
           legend=c(paste0("AUC 12:\t",AUC_12),
                    paste0("AUC 24:\t",AUC_24),
                    paste0("AUC 36:\t",AUC_36)),
           cex=0.8,
           col=colName, lwd=2, lty=1, box.lwd=1, inset=.02)
  }else{
    roc_data2=roc_data %>% set_colnames(c("TIME_TO_EVENT","EVENT","risk"))
    AUC_12=run_time_ROC(12,roc_data2,colName[1]); par(new=TRUE)
    AUC_24=run_time_ROC(24,roc_data2,colName[2]); par(new=TRUE)
    AUC_36=run_time_ROC(36,roc_data2,colName[3]); par(new=TRUE)
    legend("bottomright", 
           legend=c(paste0("AUC 12:\t",AUC_12),
                    paste0("AUC 24:\t",AUC_24),
                    paste0("AUC 36:\t",AUC_36)),
           cex=0.8,
           col=colName, lwd=2, lty=1, box.lwd=1, inset=.02)
  }
}


each_yr_roc(GES_db6_roc_data)
each_yr_roc(GES_lihc_roc_data)
each_yr_roc(GES_liri_roc_data)
each_yr_roc(IS_db6_roc_data)
each_yr_roc(IS_lihc_roc_data)
each_yr_roc(IS_liri_roc_data)


