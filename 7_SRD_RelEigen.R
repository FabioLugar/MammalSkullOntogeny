setwd("~/Onto1sp_bayesian/")
library(evolqg)
library(psych)
library(vcvComp)
library(ape)
library(plyr)
library(magrittr)
library(Matrix)
library(matrixcalc)
library(abind)
library(doParallel)
library(reshape2)
library(ggplot2)
library(cowplot)
theme_set(theme_bw())
registerDoParallel(cores=50)
matricesDist<-readRDS("matricesDist.Rds")
matricesL<-readRDS("matricesL.Rds")
info<-readRDS("ids.Rds")

matd20<-matricesL[info$gen=="Monodelphis (B)"]$d20
matd40<-matricesL[info$gen=="Monodelphis (B)"]$d40
matAdult<-matricesL[info$gen=="Monodelphis (B)"]$adults

didelphisloadings.plot<-
  rbind(data.frame(PC1=eigen(matricesL[info$gen=="Didelphis"]$zero)$vectors[,1]*-1,
                 PC2=eigen(matricesL[info$gen=="Didelphis"]$zero)$vectors[,2],
                 age="zero"),
      data.frame(PC1=eigen(matricesL[info$gen=="Didelphis"]$adults)$vectors[,1]*-1,
                 PC2=eigen(matricesL[info$gen=="Didelphis"]$adults)$vectors[,2],
                 age="adult")) %>%
  ggplot(., aes(PC1, as.numeric(age)))+
  geom_curve(aes(x=0,y=as.numeric(age), xend=PC1, yend=as.numeric(age)), 
             arrow = arrow(length = unit(0.01, "npc")),
             curvature = 0.5, alpha=0.5)+
  geom_point(shape=21)+
  geom_vline(xintercept=0)+
  scale_x_continuous(name="PC1 loadings",limits = c(-0.5,0.5))+
  scale_y_continuous(name="Age class",breaks=c(1,2),limits = c(0.5,2.5),labels = c("Zero", "Adult"))

monodelphisloadings.plot<-rbind(data.frame(PC1=eigen(matricesL[info$gen=="Monodelphis (B)"]$d20)$vectors[,1]*-1,
                 PC2=eigen(matricesL[info$gen=="Monodelphis (B)"]$d20)$vectors[,2],
                 age="d20"),
                 data.frame(PC1=eigen(matricesL[info$gen=="Monodelphis (B)"]$d40)$vectors[,1]*-1,
                            PC2=eigen(matricesL[info$gen=="Monodelphis (B)"]$d40)$vectors[,2],
                            age="d40"),
      data.frame(PC1=eigen(matricesL[info$gen=="Monodelphis (B)"]$adults)$vectors[,1]*-1,
                 PC2=eigen(matricesL[info$gen=="Monodelphis (B)"]$adults)$vectors[,2],
                 age="adult")) %>%
  ggplot(., aes(PC1, as.numeric(age)))+
  geom_curve(aes(x=0,y=as.numeric(age), xend=PC1, yend=as.numeric(age)), 
             arrow = arrow(length = unit(0.01, "npc")),
             curvature = .5, alpha=0.5)+
  geom_point(shape=21)+
  geom_vline(xintercept=0)+
  scale_x_continuous(name="PC1 loadings",limits = c(-1,0.5))+
  scale_y_continuous(name="Age class",breaks=c(1:3),limits = c(0,3.5),labels = c("d20","d40","Adult"))
                  
ggsave(plot=plot_grid(monodelphisloadings.plot,
                 didelphisloadings.plot,nrow = 2,
                 align = "hv", labels = c("A", "B")), filename = "PC1scores.pdf",width = 5,height = 8)

data.frame(Didelphis_zero=eigen(matricesL[info$gen=="Didelphis"]$zero)$vectors[,1]*-1,
           Didelphis_adult=eigen(matricesL[info$gen=="Didelphis"]$adults)$vectors[,1]*-1,
           Monodelphis_d20=eigen(matricesL[info$gen=="Monodelphis (B)"]$d20)$vectors[,1]*-1,
           Monodelphis_d40=eigen(matricesL[info$gen=="Monodelphis (B)"]$d40)$vectors[,1]*-1,
           Monodelphis_Adult=eigen(matricesL[info$gen=="Monodelphis (B)"]$adults)$vectors[,1]*-1) %>% write.csv(., file="PC1scores.csv")
          


prop.vcv.test(c(21,31),matd20,matAdult)
reld20<-relative.eigen(matd20,matAdult,method = 0)
eigen.test(c(21,31),reld20$relValues)
prop.vcv.test(c(21,31),matd40,matAdult)
reld40<-relative.eigen(matd40,matAdult,method = 0)
eigen.test(c(17,31),reld20$relValues)

pval_eigen.plot<-data.frame(comp=factor(paste(paste0("PC",1:34), 2:35, sep = "-"),levels=paste(paste0("PC",1:34), 2:35, sep = "-")),
           d20=eigen.test(c(17,31),reld20$relValues),
           d40=eigen.test(c(17,31),reld40$relValues)) %>%
  melt %>%
  ggplot(., aes(comp,value))+
  geom_line(aes(group=variable,color=variable))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=8))+
  geom_hline(yintercept = 0.05, linetype=2)+
  xlab("")+
  ylab("p-value")+
  scale_color_manual(values=c("black", "gray"))+
  annotate("label",x = 30,y = 0.05, label=expression(paste(alpha, "=0.05")),parse=TRUE)+
  ylim(0,1)
ggsave("pval_eigen.pdf",pval_eigen.plot,width = 7,height = 5) 
  
SRDd20<-SRD(matAdult,matd20)
SRDd40<-SRD(matAdult,matd40)
SRDd400<-SRD(matAdult,matricesL[info$gen=="Monodelphis (B)"]$d400)
SRD.plot<-
  rbind(data.frame(age="d20",trait=factor(rownames(SRDd20$output),rownames(SRDd20$output)),SRDd20$output),
        data.frame(age="d40",trait=factor(rownames(SRDd40$output),rownames(SRDd40$output)),SRDd40$output),
        data.frame(age="d400",trait=factor(rownames(SRDd400$output),rownames(SRDd400$output)),SRDd400$output)) %>%
  ggplot(.,
       aes(trait, ARC))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=8))+
  geom_point()+
  geom_errorbar(aes(ymin=IC.,ymax=IC..1))+
  geom_hline(aes(yintercept = ARC), 
             data.frame(age=c("d20","d40","d400"), 
                        ARC=c(median(SRDd20$output[,1]),median(SRDd40$output[,1]),median(SRDd400$output[,1]))), linetype=2)+
  ylim(c(-1,1))+
  facet_grid(age~.)+
  ylab("SRD")

ggsave("SRD.pdf",SRD.plot,width = 7,height = 6)






