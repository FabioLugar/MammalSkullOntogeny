library(evolqg)
library(plyr)
library(dplyr)
library(magrittr)
library(doParallel)
library(reshape2)
library(psych)
library(vcvComp)
library(ape)
library(Matrix)
library(matrixcalc)
library(abind)
library(ggplot2)
library(ggforce)
library(ggrepel)
library(concaveman)
library(png)
library(scales)
library(grid)
library(RCurl)
library(ggbeeswarm)
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

# Matrix repetability ----------
matricesDist<-c(llply(models$Monodelphis_D, BayesianCalculateMatrix, samples=100),
                llply(models$Monodelphis_B, BayesianCalculateMatrix, samples=100),
                llply(models$Didelphis, BayesianCalculateMatrix, samples=100),
                llply(models$Cebus, BayesianCalculateMatrix, samples=100),
                llply(models$Calomys, BayesianCalculateMatrix, samples=100))

matricesL<-llply(matricesDist, function(x) x$P)

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
    G_sagu = 0.75)

Gs<-list(G_mono=G_mono,G_calo=G_calo,G_sagu=G_sagu)

# Random Skewers -----

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

# PCOA ------

matricesLs<-matricesDist %>%
  llply(., function(x) {
    x<-x$P
    x/psych::tr(x)
  })
matricesA<-abind(matricesLs, along = 3)[,,info$age!="adults"]
meanP<-apply(matricesA,1:2,weighted.mean,w=info$ns[info$age!="adults"])
eigenP<-eigen(meanP)

dims<-20
matricesLs<-matricesDist %>%
  llply(., function(x) {
    for(i in 1:100) {
      X<-t(eigenP$vectors[,1:dims]) %*% x$Ps[i,,] %*% eigenP$vectors[,1:dims]
      x$Ps[i,1:dims,1:dims]<-X/psych::tr(X)
    }
    alply(x$Ps[1:100,1:dims,1:dims], 1, extract)
    # alply(x$Ps[1:100,,], 1, extract)
  })
matricesLs<-unlist(matricesLs[info$age!="adults"], recursive = FALSE)
eigen.phen<-MatrixDistance(matricesLs, distance = c("RiemannDist"),parallel = TRUE)

prcoa <- pr.coord(as.matrix(as.dist(eigen.phen)))

# Plots -----
theme_set(theme_bw())
shade.files <-
  dir(path = 'silhouettes/', pattern = 'shade.png', full.names = T) %>% rev(.)
colsilh<-brewer_pal(type = "qual",palette = "Set2")(5)
shade.content <- list()
for (i in 1:length(shade.files)){
  shade.content [[i]] <- rasterGrob (readPNG(shade.files [i]), interpolate = TRUE)
  # shade.content [[i]]$raster<-sub("#000000FF", colsilh[i],shade.content [[i]]$raster)
}

#---------
RS_P$corr[RS_P$corr>1]<-1

RS.plot<-
  rbind(RS_G, RS_P) %>% mutate(., Gen=factor(Gen,levels = unique(RS_G$Gen)[c(2:5,1)])) %>%
  subset(., .id!="adults") %>%
  ggplot(aes(Gen,corr))+
  geom_point(aes(size=n, shape=target_type), 
             position=position_jitterdodge(jitter.width = 0.25), alpha=1)+
  geom_vline(xintercept = c(1,3,5),
             colour = 'grey70', size = 28, alpha = 0.2) +
  theme_minimal()+
  # scale_color_brewer(type = "qual",palette = "Set2")+
  ylim(c(0.1,1))+
  coord_flip()+
  xlab("")+
  ylab("Adjusted RS similarity")+
  guides(color="none")+
  theme(panel.grid.major.y = element_blank (),
        panel.grid.major.x = element_line (color = 'black', linetype = 'dotted'),
        axis.ticks.y = element_blank ())+
  annotation_custom (shade.content [[1]], xmin = -1.75, ymin = 0.1, ymax = 0.3)+
  annotation_custom (shade.content [[2]], xmin = -3.75, ymin = 0.1, ymax = 0.3)+
  annotation_custom (shade.content [[3]], xmin = 0.1, ymin = 0.1, ymax = 0.25) +
  annotation_custom (shade.content [[4]], xmin = 2.2, ymin = 0.05, ymax = 0.25)+ 
  annotation_custom (shade.content [[5]], xmin = 4.3, ymin = 0.1, ymax = 0.3)+
  scale_size(name = "sample size")+
  scale_shape_manual(values=c(19,21), name="Target Matrix",labels=c("G","Adult P"))
RS.plot

ggsave("RS.pdf",RS.plot,width = 7,height = 5)

#----

annotation_custom2<- 
  function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data){ 
    layer(data = data, stat = StatIdentity, position = PositionIdentity, 
          geom = ggplot2:::GeomCustomAnn,
          inherit.aes = TRUE, params = list(grob = grob, 
                                            xmin = xmin, xmax = xmax, 
                                            ymin = ymin, ymax = ymax))
  }


krzsubspace.plot<-
  ggplot(melt(obsKrz), aes(as.numeric(variable), value))+
  facet_wrap(.~.id,scales = "free_y",labeller = labeller(.id = c("Monodelphis_B" = "B. Monodelphis (B)","Monodelphis_D" = "A. Monodelphis (D)","Calomys" = "E. Calomys","Cebus"= "D. Sapajus","Didelphis"="C. Didelphis")))+
  # stat_summary(aes(color=.id),fun = median, fun.min = function(x) quantile(x,0.025), fun.max = function(x) quantile(x,0.975))+
  stat_summary(fun = median, fun.min = function(x) quantile(x,0.025), fun.max = function(x) quantile(x,0.975))+
  scale_x_continuous(limits = c(0,11))+
  geom_line(aes(Var2,value, group=Var1),randKrz, linetype=2)+
  ylab("eigenvalues")+
  xlab("rank")+
  # scale_color_brewer(name="",type = "qual",palette = "Set2")+
  theme(legend.position = "none")+
  annotation_custom2(shade.content [[1]], xmin=6, xmax=11, ymin=4, data=data.frame(.id="Monodelphis_D", variable=1,value=0))+
  annotation_custom2(shade.content [[2]], xmin=6, xmax=11, ymin=3, data=data.frame(.id="Monodelphis_B", variable=1,value=0))+
  annotation_custom2(shade.content [[3]], xmin=7, xmax=11, ymin=4, data=data.frame(.id="Didelphis", variable=1,value=0))+
  annotation_custom2(shade.content [[4]], xmin=5, xmax=10, ymin=2, data=data.frame(.id="Cebus", variable=1,value=0))+
  annotation_custom2(shade.content [[5]], xmin=6, xmax=11, ymin=3.5, data=data.frame(.id="Calomys", variable=1,value=0))
krzsubspace.plot
ggsave("krzsubspace.pdf",krzsubspace.plot,width = 7,height = 5)
#----

infO<-subset(ids, age!="adults")
infO<-infO[rep(1:dim(infO)[1], each=100),]
infO$gen<-sub("Cebus","Sapajus",infO$gen)
infO$gen<-factor(infO$gen,rev(unique(infO$gen)))
means<-
  data.frame(infO,prcoa$PCoords[,]) %>% 
  group_by(gen,age) %>% summarise_all(., mean)
pvars<-var(means[,-c(1:3)])/sum(var(means[,-c(1:3)]))

pcoa.plot<-
  means %>%
  ggplot(., aes(PCo1,PCo2))+
  geom_point(aes(PCo1,PCo2,shape=gen), size=2,color="gray", data.frame(infO,prcoa$PCoords[,1:2]), show.legend = T)+
  # geom_path(aes(group=gen),means, alpha=0.5)+
  geom_text_repel(aes(label=age),means, show.legend = F,max.overlaps = 20)+
  geom_point(aes(shape=gen, fill=gen), show.legend = T)+
  coord_fixed()+
  # xlim(c(-0.5,0.5))+
  # ylim(c(-0.55,0.3))+
  scale_shape_manual(name="Ontogenetic Series",values=c(24,25,23,21,21))+
  scale_size(name="Sample size")+
  scale_fill_manual(name="Ontogenetic Series",values=c("gray","black","black","gray","black"))+
  annotation_custom (shade.content [[1]], xmin = -4, ymin = -0.8, ymax = -0.55)+   #Monodelphis
  annotation_custom (shade.content [[3]], xmin = 0.3, ymin = -0.3, ymax = 0) + #Didelphis
  annotation_custom (shade.content [[4]], xmin = 1, ymin = 0.4, ymax = 0.75)+     #Sapajus
  annotation_custom (shade.content [[5]], xmin = -6, ymin = 0.65, ymax = 0.85)+    #Calomys
  xlab(paste0("PCo1 (",round(100*prcoa$Variance$exVar[1],2),"%)"))+
  ylab(paste0("PCo2 (",round(100*prcoa$Variance$exVar[2],2),"%)"))
pcoa.plot
ggsave("pcoa.pdf",pcoa.plot,width = 8,height = 4)

pcoa34.plot<-means %>%
  ggplot(., aes(PCo3,PCo4))+
  geom_point(aes(PCo3,PCo4,shape=gen), size=2,color="gray", data.frame(infO,prcoa$PCoords), show.legend = T)+
  # geom_path(aes(group=gen),means, alpha=0.5)+
  geom_text_repel(aes(label=age),means, show.legend = F,max.overlaps = 20)+
  geom_point(aes(shape=gen, fill=gen), show.legend = T)+
  coord_fixed()+
  # xlim(c(-0.5,0.5))+
  # ylim(c(-0.55,0.3))+
  scale_shape_manual(name="Ontogenetic Series",values=c(24,25,23,21,21))+
  scale_size(name="Sample size")+
  scale_fill_manual(name="Ontogenetic Series",values=c("gray","black","black","gray","black"))+
  xlab(paste0("PCo3 (",round(100*prcoa$Variance$exVar[3],2),"%)"))+
  ylab(paste0("PCo4 (",round(100*prcoa$Variance$exVar[4],2),"%)"))
ggsave("pcoa34.pdf",pcoa34.plot,width = 7,height = 7)


#----
models<-readRDS("models.Rds")

rates<-
  ldply(1:1000, function(i) {
    models %>% ldply(., function (x) {
      sizes<-ldply(x[-length(x)], function(y)  {
        gm<-rowMeans(log(y$model$temp.Y))
        exp(rnorm(1, mean(gm), sd(gm)/sqrt(length(gm))))
      })
      age<-paste(sizes$.id[-length(sizes$.id)],sizes$.id[-1],sep = "-")
      data.frame(age_class=age,
                 rates=sizes$V1[-1]-sizes$V1[-length(sizes$V1)]) 
    }) %>% mutate(., age_fac=paste(formatC(1:26, width = 2, format = "d", flag = "0"), age_class,sep="."))
  })
rates$.id<-sub("Cebus","Sapajus",rates$.id)
rates$.id<-factor(rates$.id, levels=c("Monodelphis_D","Monodelphis_B","Didelphis","Sapajus","Calomys"))
rates$age_fac<-factor(rates$age_fac,levels = unique(rates$age_fac))
rates.plot<-
  ggplot(rates, aes(as.character(age_fac),rates))+
  facet_wrap(.~.id,scales = "free_x",dir = "v")+
  geom_jitter()+
  geom_boxplot(aes(group=age_class))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  xlab("Age class")+
  ylab("rate (%)")

ggsave("rates.pdf",rates.plot,width = 5,height = 7)

#----
matricesDist<-readRDS("matricesDist.Rds")
matricesL<-readRDS("matricesL.Rds")
matricesLs<-matricesDist %>%
  llply(., function(x) {
    x<-x$P
    x/psych::tr(x)
  })
redList<-matricesLs[ids$gen=="Monodelphis (D)" & ids$age!="adults"]
RS<-RandomSkewers(redList,parallel = T)
dists<-MatrixDistance(redList,distance = "RiemannDist",parallel = T)
df1<-data.frame(RS=c(RS$correlations),dists=c(dists), gen="Monodelphis (D)") %>%
  subset(., RS>0)
redList<-matricesL[ids$gen=="Didelphis" & ids$age!="adults"]
RS<-RandomSkewers(redList,parallel = T)
dists<-MatrixDistance(redList,distance = "RiemannDist",parallel = T)
df2<-data.frame(RS=c(RS$correlations),dists=c(dists), gen="Didelphis") %>%
  subset(., RS>0)
redList<-matricesL[ids$gen=="Cebus" & ids$age!="adults"]
RS<-RandomSkewers(redList,parallel = T)
dists<-MatrixDistance(redList,distance = "RiemannDist",parallel = T)
df3<-data.frame(RS=c(RS$correlations),dists=c(dists), gen="Cebus") %>%
  subset(., RS>0)
redList<-matricesLs[ids$gen=="Monodelphis (B)" & ids$age!="adults"]
RS<-RandomSkewers(redList,parallel = T)
dists<-MatrixDistance(redList,distance = "RiemannDist",parallel = T)
df4<-data.frame(RS=c(RS$correlations),dists=c(dists), gen="Monodelphis (B)") %>%
  subset(., RS>0)
redList<-matricesL[ids$gen=="Calomys" & ids$age!="adults"]
RS<-RandomSkewers(redList,parallel = T)
dists<-MatrixDistance(redList,distance = "RiemannDist",parallel = T)
df5<-data.frame(RS=c(RS$correlations),dists=c(dists), gen="Calomys") %>%
  subset(., RS>0)
df3$gen<-sub("Cebus","Sapajus",df3$gen)
riemRs.plot<-
  rbind(df5,df3,df2,df1,df4) %>%
  ggplot(., aes(dists, psych::fisherz(RS)))+
  geom_smooth(method="lm")+
  geom_point(aes(shape=gen, fill=gen))+
  scale_x_log10()+
  scale_y_continuous(breaks=fisherz(c(0.5,0.75, 0.875, 0.95)),labels=c(0.50,0.750, 0.875,0.950))+
  xlab("Riemannian distance")+
  ylab("Random Skewers")+
  scale_shape_manual(name="Ontogenetic Series",values=c(24,25,23,21,21))+
  scale_size(name="Sample size")+
  scale_fill_manual(name="Ontogenetic Series",values=c("gray","black","black","gray","black"))

ggsave("riemRs.pdf",riemRs.plot,width = 5,height = 4)