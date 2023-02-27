setwd("~/Onto1sp_bayesian/")
library(evolqg)
library(Matrix)
library(plyr)
library(dplyr)
library(mvtnorm)
library(doParallel)
registerDoParallel(cores=50)

models<-readRDS("models.Rds")

matricesDist<-c(llply(models$Monodelphis_D, BayesianCalculateMatrix, samples=100),
                llply(models$Monodelphis_B, BayesianCalculateMatrix, samples=100),
                llply(models$Didelphis, BayesianCalculateMatrix, samples=100),
                llply(models$Cebus, BayesianCalculateMatrix, samples=100),
                llply(models$Calomys, BayesianCalculateMatrix, samples=100))
saveRDS(matricesDist,"matricesDist.Rds")

matricesL<-llply(matricesDist, function(x) x$P)

saveRDS(matricesL,"matricesL.Rds")

modelsunlisted<-unlist(models,recursive = F)
repRS<- llply(1:length(matricesL), function(i){
  lm_f<-modelsunlisted[[i]]
  P<-matricesL[[i]]
  matrices<-llply(1:1000, function(j){
    lm_f$residuals<-rmvnorm(nrow(lm_f$residuals),sigma=P)
    BayesianCalculateMatrix(lm_f, samples=100)$P
  },.parallel = T)
  rs<-RandomSkewers(matrices, P, parallel = TRUE)
  mean(rs$correlation)
})
repRS<-unlist(repRS)
saveRDS(repRS,"repRS.RDS")

G_mono <- read.table(file = "G_mono.txt") %>% as.matrix
G_mono <- forceSymmetric(G_mono) %>% as.matrix
load('PG.Calomys.RData')
G_calo <- Calomys.Pacote $ G
G_calo <- forceSymmetric(G_calo) %>% as.matrix
G_sagu <- read.table("g-saguinus.csv") %>% as.matrix
colnames(G_sagu)<- rownames(G_sagu) <-read.table("tracos_NWM.txt")$V1
G_sagu <- G_sagu[colnames(G_calo),colnames(G_calo)]

repRSG <-
  c(G_calo = MonteCarloRep(G_calo, 20, RandomSkewers, iterations = 1000, parallel = T),
    G_mono = MonteCarloRep(G_mono, 15, RandomSkewers, iterations = 1000, parallel = T),
    G_sagu = MonteCarloRep(G_sagu, 35, RandomSkewers, iterations = 1000, parallel = T))
saveRDS(repRSG,"repRSG.RDS")
saveRDS(list(G_mono=G_mono,G_calo=G_calo,G_sagu=G_sagu),"Gs.RDS")