options(stringsAsFactors = F)
library(plyr)
library(tidyverse)
library(doParallel)
library(magrittr)
library(data.table)
library(parallel)
library(ggplot2)
library(Hmisc)



plot_nomogram=function(mydata){
  library(Hmisc)
  library(grid)
  library(lattice)
  library(Formula)
  library(ggplot2)
  library(rms)
  library(survival)



  dd=datadist(mydata)
  options(datadist="dd")
  
  formula_str=paste0("Surv(TIME_TO_EVENT,EVENT)~GES_group+IS_group+TNM")
  myformula=as.formula(formula_str)
  
  coxm <-cph(Surv(TIME_TO_EVENT,EVENT)~GES_group+IS_group+TNM,x=T,y=T,data=mydata,surv=T)
  surv<- Survival(coxm) 
  surv2<- function(x)surv(24,lp=x) 
  surv3<- function(x)surv(36,lp=x)
  plot(nomogram(coxm,fun=list(surv2,surv3),lp=F,
                funlabel=c('2-YearOS','3-YearOS'),
                maxscale=100,
                est.all=F,
                fun.at=c('0.95','0.85','0.80','0.70','0.6','0.5','0.4','0.3','0.2','0.1')),xfrac=.3)

  f<-coxph(myformula,
           x=T,y=T,data=mydata)
  sum.surv<-summary(f)
  c_index<-sum.surv$concordance
  return(list(c_index=c_index,mymodel=coxm))
}




plot_nomogram_cali_curve=function(my_formula,mydata,sampling_num,months,mycolor){
  coxm <-cph(as.formula(my_formula),
             x=T,y=T,data=mydata,surv=T)
  cal<- calibrate(coxm, cmethod='KM', method='boot', u=months, m=sampling_num,B=NROW(mydata),conf.int=T)
  plot(cal,lwd=2,lty=1,
       errbar.col=c(rgb(0,118,192,maxColorValue=255)),
       xlim=c(0,1),
	   ylim=c(0,1),
       xlab=paste0("Nomogram-PredictedProbabilityof ",months," months OS"),
       ylab=paste0("Actual ",months," months OS (proportion)"),
       col=c(rgb(192,98,83,maxColorValue=255)),
	   plot=F)
  lines(cal[,c("mean.predicted","KM")],type="b",lwd=2,col=mycolor, pch=16)
  abline(0,1,lty=3,lwd=2,col=c(rgb(0,118,192,maxColorValue=255)))
}