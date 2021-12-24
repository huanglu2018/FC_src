options(stringsAsFactors = F)
library(plyr)
library(tidyverse)
library(doParallel)
library(magrittr)
library(data.table)
library(parallel)
library(ggplot2)
library(Hmisc)
library(forestplot)

each_par=read.csv(each_csv,header = F) %>% parse_forest_info2()
outpdf_file=gsub("csv$","forest.pdf",each_csv)
pdf(file=outpdf_file,height = 5,width = 7, onefile = FALSE)
forestplot(each_par$TEXT,each_par$PARA,new_page = TRUE,
           is.summary=c(TRUE,rep(FALSE,NROW(each_par$TEXT)-1)),
           clip=c(0.05,8.0), 
           xlog=TRUE,
           align=c("l","c","c","c","c"),
           grid = structure(c(1), 
                            gp = gpar(lwd=3,lty = 1, col = "grey85")),
           hrzl_lines = list("1" = gpar(lwd=2, columns=1:4, col = "black"), 
                             "2" = gpar(lwd=2, columns=2:3, col = "black"),
                             "13" = gpar(lwd=2, columns=1:3, col = "black")),
           graph.pos=4,
           colgap=unit(0.03,"npc"),
           ci.vertices=T,
           lwd.ci=3,
           xticks=c(0.125,0.25,0.5,1,2,4,8),
           boxsize=0.2,
           graphwidth= unit(35,"mm"),
           col=fpColors(box="royalblue",line="red4", summary="royalblue")
dev.off()


each_par=read.csv(each_csv,header = F) %>% parse_forest_info2()
outpdf_file=gsub("csv$","forest.pdf",each_csv)
pdf(file=outpdf_file,height = 5,width = 7, onefile = FALSE)
forestplot(each_par$TEXT,each_par$PARA,new_page = TRUE,
           is.summary=c(TRUE,rep(FALSE,NROW(each_par$TEXT)-1)),
           clip=c(0.04,8.0), 
           xlog=TRUE,
           align=c("l","c","c","c","c"),
           grid = structure(c(1), 
                            gp = gpar(lwd=3,lty = 1, col = "grey85")),
           hrzl_lines = list("1" = gpar(lwd=2, columns=1:4, col = "black"), 
                             "2" = gpar(lwd=2, columns=2:3, col = "black"),
                             "13" = gpar(lwd=2, columns=1:3, col = "black")),
           graph.pos=4,
           colgap=unit(0.03,"npc"),
           ci.vertices=T,
           lwd.ci=3,
           xticks=c(0.125,0.25,0.5,1,2,4,8),
           boxsize=0.2,
           graphwidth= unit(35,"mm"),
           col=fpColors(box="royalblue",line="red4", summary="royalblue")
dev.off()
