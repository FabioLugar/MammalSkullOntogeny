setwd("~/Onto1sp_bayesian/")
library(evolqg)
library(Matrix)
library(plyr)
library(dplyr)
library(mvtnorm)
library(doParallel)
registerDoParallel(cores=50)

models<-readRDS("models.Rds")
Gs<-readRDS("Gs.RDS")
matricesL<-readRDS("matricesL.Rds")
matricesDist<-readRDS("matricesDist.Rds")

ns<-c(llply(models$Monodelphis_D, function(x){
  nrow(x$residuals)}),
  llply(models$Monodelphis_B, function(x){
    nrow(x$residuals)}),
  llply(models$Didelphis,  function(x){
    nrow(x$residuals)}),
  llply(models$Cebus,  function(x){
    nrow(x$residuals)}),
  llply(models$Calomys, function(x){
    nrow(x$residuals)}))
ns<-unlist(ns)

repRS<-readRDS("repRS.RDS")
repRSG<-readRDS("repRSG.RDS")
repRSG[3]<-0.75
age<-names(matricesL)
gen<-c(rep("Monodelphis (D)",length(models$Monodelphis_D)),
       rep("Monodelphis (B)",length(models$Monodelphis_B)),
       rep("Didelphis",length(models$Didelphis)),
       rep("Cebus",length(models$Cebus)),
       rep("Calomys",length(models$Calomys)))

RS_G<-rbind(
  data.frame(Gen="Calomys", age_type="B", RandomSkewers(matricesL[gen=="Calomys"],Gs$G_calo)) %>% 
    dplyr::select(.,Gen,.id,age_type,correlation) %>% rename(., raw=correlation) %>%  
    mutate(., corr=raw/sqrt(repRSG["G_calo"]*repRS[gen=="Calomys"]),
           n=ns[gen=="Calomys"],
           target_type="G"),
  data.frame(Gen="Monodelphis (D)", age_type="D", RandomSkewers(matricesL[gen=="Monodelphis (D)"],Gs$G_mono)) %>% 
    dplyr::select(.,Gen,.id,age_type,correlation) %>% rename(., raw=correlation) %>%  
    mutate(., corr=raw/sqrt(repRSG["G_mono"]*repRS[gen=="Monodelphis (D)"]),
           n=ns[gen=="Monodelphis (D)"],
           target_type="G"),
  data.frame(Gen="Monodelphis (B)", age_type="B", RandomSkewers(matricesL[gen=="Monodelphis (B)"],Gs$G_mono)) %>% 
    dplyr::select(.,Gen,.id,age_type,correlation) %>% rename(., raw=correlation) %>%  
    mutate(., corr=raw/sqrt(repRSG["G_mono"]*repRS[gen=="Monodelphis (B)"]),
           n=ns[gen=="Monodelphis (B)"],
           target_type="G"),
  data.frame(Gen="Didelphis", age_type="D", RandomSkewers(matricesL[gen=="Didelphis"],Gs$G_mono)) %>% 
    dplyr::select(.,Gen,.id,age_type,correlation) %>% rename(., raw=correlation) %>%  
    mutate(., corr=raw/sqrt(repRSG["G_mono"]*repRS[gen=="Didelphis"]),
           n=ns[gen=="Didelphis"],
           target_type="G"),
  data.frame(Gen="Sapajus", age_type="D", RandomSkewers(matricesL[gen=="Cebus"],Gs$G_sagu)) %>% 
    dplyr::select(.,Gen,.id,age_type,correlation) %>% rename(., raw=correlation) %>%  
    mutate(., corr=raw/sqrt(repRSG["G_sagu"]*repRS[gen=="Cebus"]),
           n=ns[gen=="Cebus"],
           target_type="G"))

RS_P<-rbind(
  data.frame(Gen="Calomys", age_type="B", RandomSkewers(matricesL[gen=="Calomys"][-sum(gen=="Calomys")],matricesL[gen=="Calomys"][["adults"]])) %>% 
    dplyr::select(.,Gen,.id,age_type,correlation) %>% rename(., raw=correlation) %>%  
    mutate(., corr=raw/sqrt(repRS[gen=="Calomys"][sum(gen=="Calomys")]*repRS[gen=="Calomys"][-sum(gen=="Calomys")]),
           n=ns[gen=="Calomys"][-sum(gen=="Calomys")],
           target_type="P"),
  data.frame(Gen="Monodelphis (D)", age_type="B", RandomSkewers(matricesL[gen=="Monodelphis (D)"][-sum(gen=="Monodelphis (D)")],matricesL[gen=="Monodelphis (D)"][["adults"]])) %>% 
    dplyr::select(.,Gen,.id,age_type,correlation) %>% rename(., raw=correlation) %>%  
    mutate(., corr=raw/sqrt(repRS[gen=="Monodelphis (D)"][sum(gen=="Monodelphis (D)")]*repRS[gen=="Monodelphis (D)"][-sum(gen=="Monodelphis (D)")]),
           n=ns[gen=="Monodelphis (D)"][-sum(gen=="Monodelphis (D)")],
           target_type="P"),
  data.frame(Gen="Monodelphis (B)", age_type="B", RandomSkewers(matricesL[gen=="Monodelphis (B)"][-sum(gen=="Monodelphis (B)")],matricesL[gen=="Monodelphis (B)"][["adults"]])) %>% 
    dplyr::select(.,Gen,.id,age_type,correlation) %>% rename(., raw=correlation) %>%  
    mutate(., corr=raw/sqrt(repRS[gen=="Monodelphis (B)"][sum(gen=="Monodelphis (B)")]*repRS[gen=="Monodelphis (B)"][-sum(gen=="Monodelphis (B)")]),
           n=ns[gen=="Monodelphis (B)"][-sum(gen=="Monodelphis (B)")],
           target_type="P"),
  data.frame(Gen="Didelphis", age_type="B", RandomSkewers(matricesL[gen=="Didelphis"][-sum(gen=="Didelphis")],matricesL[gen=="Didelphis"][["adults"]])) %>% 
    dplyr::select(.,Gen,.id,age_type,correlation) %>% rename(., raw=correlation) %>%  
    mutate(., corr=raw/sqrt(repRS[gen=="Didelphis"][sum(gen=="Didelphis")]*repRS[gen=="Didelphis"][-sum(gen=="Didelphis")]),
           n=ns[gen=="Didelphis"][-sum(gen=="Didelphis")],
           target_type="P"),
  data.frame(Gen="Sapajus", age_type="B", RandomSkewers(matricesL[gen=="Cebus"][-sum(gen=="Cebus")],matricesL[gen=="Cebus"][["adults"]])) %>% 
    dplyr::select(.,Gen,.id,age_type,correlation) %>% rename(., raw=correlation) %>%  
    mutate(., corr=raw/sqrt(repRS[gen=="Cebus"][sum(gen=="Cebus")]*repRS[gen=="Cebus"][-sum(gen=="Cebus")]),
           n=ns[gen=="Cebus"][-sum(gen=="Cebus")],
           target_type="P"))

saveRDS(RS_P,"RS_P.Rds")
saveRDS(RS_G,"RS_G.Rds")
saveRDS(data.frame(gen,age,ns),"ids.Rds")
