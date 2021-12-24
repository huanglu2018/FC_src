options(stringsAsFactors = F)
library(plyr)
library(tidyverse)
library(doParallel)
library(magrittr)
library(data.table)
library(parallel)
library(ggplot2)
library(Hmisc)



limma_de_num=function(mytest,mydesign,adjPvalueCutoff,FC_top_num,mycoef){
  library(limma)
  fit <- lmFit(mytest, mydesign)
  fit <- eBayes(fit)
  allGeneSets <- topTable(fit, coef=mycoef, number=Inf)
  
  pfiltGeneSets <- topTable(fit, coef=mycoef,p.value=adjPvalueCutoff, number=Inf,adjust="fdr")
  
  FC_threshold=ifelse(NROW(pfiltGeneSets) < FC_top_num,
                      sort(abs(pfiltGeneSets$logFC),decreasing = T)[NROW(pfiltGeneSets)]-1e-9,
                      sort(abs(pfiltGeneSets$logFC),decreasing = T)[FC_top_num]-1e-9)
  print(paste0("FC threshold: ",FC_threshold))
  DEgeneSets <- topTable(fit, coef=mycoef, number=Inf,
                         p.value=adjPvalueCutoff,
                         lfc=FC_threshold,
                         adjust="fdr")
  return(list(all=allGeneSets,DE=DEgeneSets))
}

res_limma=limma_de_num(mytest,mydesign,adjPvalueCutoff,40,colname2)
T_ALLres=res_limma$all
T_DEres_raw=res_limma$DE


T_DEres=T_DEres_raw %>%
  rownames_to_column(var="terms") %>%
  mutate(term_simp=gsub("^GO_","",terms) %>% 
  		gsub("^KEGG_","",.) %>% 
		gsub("_PROCESS$","",.) %>% 
		gsub("_PATHWAY$","",.)) %>%
  mutate(term_simp=tolower(term_simp)) %>%
  mutate(code=1:NROW(T_DEres_raw)) %>%
  mutate(term_simp_code=paste0(term_simp," (",code,") ")) %>%
  column_to_rownames(var="terms")

for (i in rownames(T_DEres)) {
  rownames(T_ALLres)[match(i,rownames(T_ALLres))]=T_DEres[match(i,rownames(T_DEres)),"code"]
}

# heatmap
heat_data=mytest %>% as.data.frame() %>% 
  .[rownames(T_DEres_raw),IS_order_ind]%>% 
  rownames_to_column(var="terms") %>% 
  mutate(term_simp=tolower(term_simp)) %>% 
  mutate(code=1:NROW(T_DEres_raw)) %>% 
  mutate(term_simp_code=paste0(term_simp," (",code,") ")) %>% 
  mutate(term_simp_code=gsub("_"," ",term_simp_code)) %>% 
  dplyr::select(-c(term_simp,code,terms)) %>% 
  column_to_rownames(var="term_simp_code")

t1=pheatmap::pheatmap(heat_data,
                      # color = colorRampPalette(c("steelblue", "white","tomato"))(40),
                      cluster_rows = T, 
                      cluster_cols = F,
                      clustering_method="complete",
                      fontsize_row = 8,height = 11,
                      annotation_col = risk_imm_info_cmb[,c(groupmethod),drop=F],
                      show_colnames = F)
T1 = ggplotify::as.ggplot(t1)

