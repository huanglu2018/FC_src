options(stringsAsFactors = F)
library(Hmisc)
library(grid)
library(lattice)
library(Formula)
library(ggplot2)
library(rms)
library(survival)
library(MASS)
library(survival)
library(magrittr)


# calibration curves

my_formula="Surv(TIME_TO_EVENT,EVENT) ~ GES_group + TNM + IS_group"
plot_nomogram_cali_curve(my_formula,mydata,40,24,colName[2]); par(new=TRUE)
plot_nomogram_cali_curve(my_formula,mydata,40,32,colName[3])


# deciscion curves 

nomo_score=predict(plot_nomogram(mydata)$mymodel) %>%
  as.data.frame() %>% 
  set_colnames("nomo_score")
nomodata.nomo_score=merge(mydata,nomo_score,by=0) %>% 
  column_to_rownames(var="Row.names")

nomodata.dca=nomodata.nomo_score
nomodata.dca[nomodata.dca=="risk high"]=1
nomodata.dca[nomodata.dca=="risk low"]=0
nomodata.dca[nomodata.dca=="I"]=1
nomodata.dca[nomodata.dca=="II"]=2
nomodata.dca[nomodata.dca=="III"]=3
for (i in 1:NCOL(nomodata.dca)){
  nomodata.dca[,i]=as.numeric(nomodata.dca[,i])
}

source('dca.R')

result2=stdca(data=nomodata.dca, 
              outcome="EVENT", 
              ttoutcome="TIME_TO_EVENT",
              timepoint=24, 
              predictors=c("TNM","IS_group","GES_group","nomo_score"), 
              probability=c(F,F,F,F), 
              xstop=0.75)
result3=stdca(data=nomodata.dca, 
              outcome="EVENT", 
              ttoutcome="TIME_TO_EVENT",
              timepoint=36, 
              predictors=c("TNM","IS_group","GES_group","nomo_score"), 
              probability=c(F,F,F,F), 
              xstop=1)
