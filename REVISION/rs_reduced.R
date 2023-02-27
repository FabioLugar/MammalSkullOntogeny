setwd("~/Onto1sp_bayesian/REVISION/")
library(evolqg)
library(Matrix)
library(plyr)
library(dplyr)
library(mvtnorm)
library(doParallel)
registerDoParallel(cores=50)
parallel::mcaffinity(1:50)
models<-readRDS(file = "~/Onto1sp_bayesian/models.Rds")
Gs<-readRDS("~/Onto1sp_bayesian/Gs.RDS")
matricesL<-readRDS("~/Onto1sp_bayesian/matricesL.Rds")


#fifteen-------
allTraits<-colnames(Gs$G_calo)
colnames(Gs$G_mono)<-rownames(Gs$G_mono)<-allTraits
fifteenTraits<-c("NSLNA","NSLZS","NSLZI","ISPNS","PTBA","NABR","PTEAM","PTAS","ZSZI","PMZS","PMZI","PTZYGO","ZIZYGO","BRLD","PMMT")
models15<-models
Gs15<-Gs
matricesL15<-matricesL

for(i in seq_along(matricesL)) {
  matricesL15[[i]]<- matricesL[[i]][fifteenTraits,fifteenTraits]
}

for(i in seq_along(Gs)) {
  Gs15[[i]]<- Gs[[i]][fifteenTraits,fifteenTraits]
}


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

repRS<-readRDS("~/Onto1sp_bayesian/repRS.RDS")
repRSG<-readRDS("~/Onto1sp_bayesian/repRSG.RDS")
repRSG[3]<-0.75

age<-names(matricesL)
gen<-c(rep("Monodelphis (D)",length(models$Monodelphis_D)),
       rep("Monodelphis (B)",length(models$Monodelphis_B)),
       rep("Didelphis",length(models$Didelphis)),
       rep("Cebus",length(models$Cebus)),
       rep("Calomys",length(models$Calomys)))

RS_G<-rbind(
  data.frame(Gen="Calomys", age_type="B", RandomSkewers(matricesL15[gen=="Calomys"],Gs15$G_calo)) %>% 
    dplyr::select(.,Gen,.id,age_type,correlation) %>% rename(., raw=correlation) %>%  
    mutate(., corr=raw/sqrt(repRSG["G_calo"]*repRS[gen=="Calomys"]),
           n=ns[gen=="Calomys"],
           target_type="G"),
  data.frame(Gen="Monodelphis (D)", age_type="D", RandomSkewers(matricesL15[gen=="Monodelphis (D)"],Gs15$G_mono)) %>% 
    dplyr::select(.,Gen,.id,age_type,correlation) %>% rename(., raw=correlation) %>%  
    mutate(., corr=raw/sqrt(repRSG["G_mono"]*repRS[gen=="Monodelphis (D)"]),
           n=ns[gen=="Monodelphis (D)"],
           target_type="G"),
  data.frame(Gen="Monodelphis (B)", age_type="B", RandomSkewers(matricesL15[gen=="Monodelphis (B)"],Gs15$G_mono)) %>% 
    dplyr::select(.,Gen,.id,age_type,correlation) %>% rename(., raw=correlation) %>%  
    mutate(., corr=raw/sqrt(repRSG["G_mono"]*repRS[gen=="Monodelphis (B)"]),
           n=ns[gen=="Monodelphis (B)"],
           target_type="G"),
  data.frame(Gen="Didelphis", age_type="D", RandomSkewers(matricesL15[gen=="Didelphis"],Gs15$G_mono)) %>% 
    dplyr::select(.,Gen,.id,age_type,correlation) %>% rename(., raw=correlation) %>%  
    mutate(., corr=raw/sqrt(repRSG["G_mono"]*repRS[gen=="Didelphis"]),
           n=ns[gen=="Didelphis"],
           target_type="G"),
  data.frame(Gen="Sapajus", age_type="D", RandomSkewers(matricesL15[gen=="Cebus"],Gs15$G_sagu)) %>% 
    dplyr::select(.,Gen,.id,age_type,correlation) %>% rename(., raw=correlation) %>%  
    mutate(., corr=raw/sqrt(repRSG["G_sagu"]*repRS[gen=="Cebus"]),
           n=ns[gen=="Cebus"],
           target_type="G"))

RS_P<-rbind(
  data.frame(Gen="Calomys", age_type="B", RandomSkewers(matricesL15[gen=="Calomys"][-sum(gen=="Calomys")],matricesL15[gen=="Calomys"][["adults"]])) %>% 
    dplyr::select(.,Gen,.id,age_type,correlation) %>% rename(., raw=correlation) %>%  
    mutate(., corr=raw/sqrt(repRS[gen=="Calomys"][sum(gen=="Calomys")]*repRS[gen=="Calomys"][-sum(gen=="Calomys")]),
           n=ns[gen=="Calomys"][-sum(gen=="Calomys")],
           target_type="P"),
  data.frame(Gen="Monodelphis (D)", age_type="B", RandomSkewers(matricesL15[gen=="Monodelphis (D)"][-sum(gen=="Monodelphis (D)")],matricesL15[gen=="Monodelphis (D)"][["adults"]])) %>% 
    dplyr::select(.,Gen,.id,age_type,correlation) %>% rename(., raw=correlation) %>%  
    mutate(., corr=raw/sqrt(repRS[gen=="Monodelphis (D)"][sum(gen=="Monodelphis (D)")]*repRS[gen=="Monodelphis (D)"][-sum(gen=="Monodelphis (D)")]),
           n=ns[gen=="Monodelphis (D)"][-sum(gen=="Monodelphis (D)")],
           target_type="P"),
  data.frame(Gen="Monodelphis (B)", age_type="B", RandomSkewers(matricesL15[gen=="Monodelphis (B)"][-sum(gen=="Monodelphis (B)")],matricesL15[gen=="Monodelphis (B)"][["adults"]])) %>% 
    dplyr::select(.,Gen,.id,age_type,correlation) %>% rename(., raw=correlation) %>%  
    mutate(., corr=raw/sqrt(repRS[gen=="Monodelphis (B)"][sum(gen=="Monodelphis (B)")]*repRS[gen=="Monodelphis (B)"][-sum(gen=="Monodelphis (B)")]),
           n=ns[gen=="Monodelphis (B)"][-sum(gen=="Monodelphis (B)")],
           target_type="P"),
  data.frame(Gen="Didelphis", age_type="B", RandomSkewers(matricesL15[gen=="Didelphis"][-sum(gen=="Didelphis")],matricesL15[gen=="Didelphis"][["adults"]])) %>% 
    dplyr::select(.,Gen,.id,age_type,correlation) %>% rename(., raw=correlation) %>%  
    mutate(., corr=raw/sqrt(repRS[gen=="Didelphis"][sum(gen=="Didelphis")]*repRS[gen=="Didelphis"][-sum(gen=="Didelphis")]),
           n=ns[gen=="Didelphis"][-sum(gen=="Didelphis")],
           target_type="P"),
  data.frame(Gen="Sapajus", age_type="B", RandomSkewers(matricesL15[gen=="Cebus"][-sum(gen=="Cebus")],matricesL15[gen=="Cebus"][["adults"]])) %>% 
    dplyr::select(.,Gen,.id,age_type,correlation) %>% rename(., raw=correlation) %>%  
    mutate(., corr=raw/sqrt(repRS[gen=="Cebus"][sum(gen=="Cebus")]*repRS[gen=="Cebus"][-sum(gen=="Cebus")]),
           n=ns[gen=="Cebus"][-sum(gen=="Cebus")],
           target_type="P"))

saveRDS(RS_P,"RS_P15.Rds")
saveRDS(RS_G,"RS_G15.Rds")
saveRDS(data.frame(gen,age,ns),"ids.Rds")

#five-------
fiveTraits<-c("NSLNA","ISPNS","ZIZYGO","PTAS","BRLD")
models5<-models
Gs5<-Gs
matricesL5<-matricesL

for(i in seq_along(matricesL)) {
  matricesL5[[i]]<- matricesL[[i]][fiveTraits,fiveTraits]
}

for(i in seq_along(Gs)) {
  Gs5[[i]]<- Gs[[i]][fiveTraits,fiveTraits]
}


RS_G<-rbind(
  data.frame(Gen="Calomys", age_type="B", RandomSkewers(matricesL5[gen=="Calomys"],Gs5$G_calo)) %>% 
    dplyr::select(.,Gen,.id,age_type,correlation) %>% rename(., raw=correlation) %>%  
    mutate(., corr=raw/sqrt(repRSG["G_calo"]*repRS[gen=="Calomys"]),
           n=ns[gen=="Calomys"],
           target_type="G"),
  data.frame(Gen="Monodelphis (D)", age_type="D", RandomSkewers(matricesL5[gen=="Monodelphis (D)"],Gs5$G_mono)) %>% 
    dplyr::select(.,Gen,.id,age_type,correlation) %>% rename(., raw=correlation) %>%  
    mutate(., corr=raw/sqrt(repRSG["G_mono"]*repRS[gen=="Monodelphis (D)"]),
           n=ns[gen=="Monodelphis (D)"],
           target_type="G"),
  data.frame(Gen="Monodelphis (B)", age_type="B", RandomSkewers(matricesL5[gen=="Monodelphis (B)"],Gs5$G_mono)) %>% 
    dplyr::select(.,Gen,.id,age_type,correlation) %>% rename(., raw=correlation) %>%  
    mutate(., corr=raw/sqrt(repRSG["G_mono"]*repRS[gen=="Monodelphis (B)"]),
           n=ns[gen=="Monodelphis (B)"],
           target_type="G"),
  data.frame(Gen="Didelphis", age_type="D", RandomSkewers(matricesL5[gen=="Didelphis"],Gs5$G_mono)) %>% 
    dplyr::select(.,Gen,.id,age_type,correlation) %>% rename(., raw=correlation) %>%  
    mutate(., corr=raw/sqrt(repRSG["G_mono"]*repRS[gen=="Didelphis"]),
           n=ns[gen=="Didelphis"],
           target_type="G"),
  data.frame(Gen="Sapajus", age_type="D", RandomSkewers(matricesL5[gen=="Cebus"],Gs5$G_sagu)) %>% 
    dplyr::select(.,Gen,.id,age_type,correlation) %>% rename(., raw=correlation) %>%  
    mutate(., corr=raw/sqrt(repRSG["G_sagu"]*repRS[gen=="Cebus"]),
           n=ns[gen=="Cebus"],
           target_type="G"))

RS_P<-rbind(
  data.frame(Gen="Calomys", age_type="B", RandomSkewers(matricesL5[gen=="Calomys"][-sum(gen=="Calomys")],matricesL5[gen=="Calomys"][["adults"]])) %>% 
    dplyr::select(.,Gen,.id,age_type,correlation) %>% rename(., raw=correlation) %>%  
    mutate(., corr=raw/sqrt(repRS[gen=="Calomys"][sum(gen=="Calomys")]*repRS[gen=="Calomys"][-sum(gen=="Calomys")]),
           n=ns[gen=="Calomys"][-sum(gen=="Calomys")],
           target_type="P"),
  data.frame(Gen="Monodelphis (D)", age_type="B", RandomSkewers(matricesL5[gen=="Monodelphis (D)"][-sum(gen=="Monodelphis (D)")],matricesL5[gen=="Monodelphis (D)"][["adults"]])) %>% 
    dplyr::select(.,Gen,.id,age_type,correlation) %>% rename(., raw=correlation) %>%  
    mutate(., corr=raw/sqrt(repRS[gen=="Monodelphis (D)"][sum(gen=="Monodelphis (D)")]*repRS[gen=="Monodelphis (D)"][-sum(gen=="Monodelphis (D)")]),
           n=ns[gen=="Monodelphis (D)"][-sum(gen=="Monodelphis (D)")],
           target_type="P"),
  data.frame(Gen="Monodelphis (B)", age_type="B", RandomSkewers(matricesL5[gen=="Monodelphis (B)"][-sum(gen=="Monodelphis (B)")],matricesL5[gen=="Monodelphis (B)"][["adults"]])) %>% 
    dplyr::select(.,Gen,.id,age_type,correlation) %>% rename(., raw=correlation) %>%  
    mutate(., corr=raw/sqrt(repRS[gen=="Monodelphis (B)"][sum(gen=="Monodelphis (B)")]*repRS[gen=="Monodelphis (B)"][-sum(gen=="Monodelphis (B)")]),
           n=ns[gen=="Monodelphis (B)"][-sum(gen=="Monodelphis (B)")],
           target_type="P"),
  data.frame(Gen="Didelphis", age_type="B", RandomSkewers(matricesL5[gen=="Didelphis"][-sum(gen=="Didelphis")],matricesL5[gen=="Didelphis"][["adults"]])) %>% 
    dplyr::select(.,Gen,.id,age_type,correlation) %>% rename(., raw=correlation) %>%  
    mutate(., corr=raw/sqrt(repRS[gen=="Didelphis"][sum(gen=="Didelphis")]*repRS[gen=="Didelphis"][-sum(gen=="Didelphis")]),
           n=ns[gen=="Didelphis"][-sum(gen=="Didelphis")],
           target_type="P"),
  data.frame(Gen="Sapajus", age_type="B", RandomSkewers(matricesL5[gen=="Cebus"][-sum(gen=="Cebus")],matricesL5[gen=="Cebus"][["adults"]])) %>% 
    dplyr::select(.,Gen,.id,age_type,correlation) %>% rename(., raw=correlation) %>%  
    mutate(., corr=raw/sqrt(repRS[gen=="Cebus"][sum(gen=="Cebus")]*repRS[gen=="Cebus"][-sum(gen=="Cebus")]),
           n=ns[gen=="Cebus"][-sum(gen=="Cebus")],
           target_type="P"))

saveRDS(RS_P,"RS_P5.Rds")
saveRDS(RS_G,"RS_G5.Rds")

