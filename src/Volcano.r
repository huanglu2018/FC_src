
options(stringsAsFactors = F)
library(plyr)
library(tidyverse)
library(doParallel)
library(magrittr)
library(data.table)
library(parallel)
library(ggplot2)
library(Hmisc)

draw_volcano=function(Dat,adjPvalueCutoff,logFCcutoff,myxlim){
  library(ggplot2)
  library(ggrepel)
  library(ggThemeAssist)
  Dat$threshold = factor(ifelse(Dat$adj.P.Val < adjPvalueCutoff & abs(Dat$logFC) >= logFCcutoff, 
                                ifelse(Dat$logFC>= logFCcutoff ,'Up','Down'),
                                'NoSignifi'),
                         levels=c('Up','Down','NoSignifi'))
  
  Dat2=Dat %>% rownames_to_column(var="id")
  
  p=ggplot(Dat2,aes(x=logFC,y=-log10(adj.P.Val),color=threshold))+
    geom_point()+
    scale_color_manual(values=c("#DC143C","#00008B","#808080"))+#color
    geom_text_repel(
      data = Dat2[Dat2$adj.P.Val<adjPvalueCutoff&abs(Dat2$logFC)>logFCcutoff,],
      aes(label = id),
      size = 2.5,
      segment.color = "black", show.legend = FALSE )+
    theme_bw()+
    theme(
      legend.title = element_blank(),
      panel.grid.major =element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank()
    )+
    ylab('-log10 (p-adj)')+
    xlab('log2 (FoldChange)')+
    geom_vline(xintercept=c(logFCcutoff*-1,logFCcutoff),lty=3,col="black",lwd=0.5) +
    geom_hline(yintercept = -log10(adjPvalueCutoff),lty=3,col="black",lwd=0.5) + 
    xlim(myxlim*-1,myxlim)
  return(p)
}

