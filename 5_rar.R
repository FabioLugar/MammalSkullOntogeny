setwd("~/Onto1sp_bayesian/")
library(evolqg)
library(psych)
library(vcvComp)
library(ape)
library(plyr)
library(dplyr)
library(magrittr)
library(Matrix)
library(matrixcalc)
library(abind)
library(doParallel)
library(reshape2)
library(ggplot2)
registerDoParallel(cores=50)
theme_set(theme_bw())

models<-readRDS("models.Rds")
adult.model<-models$Monodelphis_B$adults
adult.mat<-BayesianCalculateMatrix(adult.model, samples = 100) 

rarvalues<-
  ldply(seq(30, 10, by = -2), function(n){
    laply(1:100, function(i){
      lm.temp<-adult.model
      X<-sample_n(data.frame(adult.model$residuals), n)
      lm.temp<-lm(as.matrix(X)~1)
      adult.mat.rar<-BayesianCalculateMatrix(lm.temp, samples = 100) 
      RandomSkewers(adult.mat.rar$P,adult.mat$P)[1]
    },.parallel = T)
  })


empRS<-c(RandomSkewers(BayesianCalculateMatrix(models$Monodelphis_B$adults, samples = 100)$P,
              BayesianCalculateMatrix(models$Monodelphis_B$d20, samples = 100)$P)[1],
         RandomSkewers(BayesianCalculateMatrix(models$Monodelphis_B$adults, samples = 100)$P,
                       BayesianCalculateMatrix(models$Monodelphis_B$d40, samples = 100)$P)[1])


rar.plot<-
  data.frame(Sample=seq(30, 10, by = -2),rarvalues) %>% melt(., id.vars="Sample") %>%
  ggplot(., aes(Sample, value))+
  stat_summary(fun = median,
               fun.max = function(x) quantile(x,0.025),
               fun.min = function(x) quantile(x,0.975))+
  ylim(0.5,1)+
  ylab("Random Skewer")+
  geom_hline(yintercept=empRS, linetype=2)+
  annotate("label",x =dim(models$Monodelphis_B$d20$residuals)[1] ,y = empRS[1], label="d20")+
  annotate("label",x =dim(models$Monodelphis_B$d40$residuals)[1] ,y = empRS[2], label="d40")
ggsave("rar.pdf",rar.plot,width = 5,height = 5)
