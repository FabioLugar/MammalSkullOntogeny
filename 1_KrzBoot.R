setwd("~/Onto1sp_bayesian/")
library(evolqg)
library(plyr)
library(dplyr)
library(magrittr)
library(doParallel)
library(reshape2)
registerDoParallel(cores=50)
models<-readRDS(file = "models.Rds")

# KrzBoot----------

krzsubspaceboot<-
  llply(models, function(x) {
    evolqg::KrzSubspaceBootstrap(x[-length(x)], rep = 1000, MCMCsamples = 100, parallel = TRUE) # excluding the "adult" matrix
    })

obsKrz<-ldply(krzsubspaceboot,extract2, 1)
obsKrz$.id<-factor(obsKrz$.id,levels = unique(obsKrz$.id))
randKrz<-llply(krzsubspaceboot,function(x) {
  x<-x[[2]]
  apply(x,2,quantile,probs=c(0.025,0.975))
}) %>% melt
names(randKrz)[4]<-".id"
randKrz$.id<-factor(randKrz$.id,levels = unique(randKrz$.id))


saveRDS(randKrz,"randKrz.Rds")
saveRDS(obsKrz,"obsKrz.Rds")