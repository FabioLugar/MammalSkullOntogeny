setwd("~/Onto1sp_bayesian/REVISION/")
library(evolqg)
library(plyr)
library(dplyr)
library(magrittr)
library(doParallel)
library(reshape2)
registerDoParallel(cores=30)
parallel::mcaffinity(51:80)
models<-readRDS(file = "~/Onto1sp_bayesian/models.Rds")

#fifteen
fifteenTraits<-c("NSLNA","NSLZS","NSLZI","ISPNS","PTBA","NABR","PTEAM","PTAS","ZSZI","PMZS","PMZI","PTZYGO","ZIZYGO","BRLD","PMMT")
models15<-models
for(i in seq_along(models)) for(j in seq_along(models[[i]])) {
  models15[[i]][[j]]<- lm(models[[i]][[j]]$res[,fifteenTraits]~1)
}
krzsubspaceboot<-
  llply(models15, function(x) {
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


saveRDS(randKrz,"randKrz15.Rds")
saveRDS(obsKrz,"obsKrz15.Rds")

#five traits
fiveTraits<-c("NSLNA","ISPNS","ZIZYGO","PTAS","BRLD")
models5<-models
for(i in seq_along(models)) for(j in seq_along(models[[i]])) {
  models5[[i]][[j]]<- lm(models[[i]][[j]]$res[,fiveTraits]~1)
}

krzsubspaceboot<-
  llply(models5, function(x) {
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


saveRDS(randKrz,"randKrz5.Rds")
saveRDS(obsKrz,"obsKrz5.Rds")

