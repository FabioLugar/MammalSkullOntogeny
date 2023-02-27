setwd("~/Onto1sp_bayesian/REVISION/")
library(plyr)
library(evolqg)
library(dplyr)
library(magrittr)
library(ggplot2)
library(ggforce)
library(psych)
library(ggrepel)
library(reshape2)
library(concaveman)
library(png)
library(scales)
library(grid)
library(RCurl)
library(doParallel)
library(ggbeeswarm)
registerDoParallel(cores=30)
theme_set(theme_bw())
shade.files <-
  dir(path = '../silhouettes/', pattern = 'shade.png', full.names = T) %>% rev(.)
colsilh<-brewer_pal(type = "qual",palette = "Set2")(5)
shade.content <- list()
for (i in 1:length(shade.files)){
  shade.content [[i]] <- rasterGrob (readPNG(shade.files [i]), interpolate = TRUE)
  # shade.content [[i]]$raster<-sub("#000000FF", colsilh[i],shade.content [[i]]$raster)
}

annotation_custom2<-
  function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data){
    layer(data = data, stat = StatIdentity, position = PositionIdentity,
          geom = ggplot2:::GeomCustomAnn,
          inherit.aes = TRUE, params = list(grob = grob,
                                            xmin = xmin, xmax = xmax,
                                            ymin = ymin, ymax = ymax))
  }

#RS 15traits--------
RS_G<-readRDS("RS_G15.Rds")
RS_P<-readRDS("RS_P15.Rds")
RS_P$corr[RS_P$corr>1]<-1
# show_col(names(sort(table(shade.content[[1]]$raster))))

# shade.content[[3]]$raster[grep("00000000",shade.content[[3]]$raster)]<-colsilh[3]
# cols<-alply(names(sort(table(shade.content[[3]]$raster))),1, color.id)

RS.plot<-
  rbind(RS_G, RS_P) %>% mutate(., Gen=factor(Gen,levels = unique(RS_G$Gen)[c(2:5,1)])) %>%
  subset(., .id!="adults") %>%
  ggplot(aes(Gen,corr))+
  geom_point(aes(size=n, shape=target_type), 
             position=position_jitterdodge(jitter.width = 0.25), alpha=1)+
  # geom_label_repel(aes(label=.id, fill=target_type),position=position_jitterdodge(jitter.width = 0.25))+
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
  scale_shape_manual(values=c(19,21), name="Target Matrix",labels=c("G","Adult P"))+
  scale_fill_manual(values=c("gray","white"), name="Target Matrix",labels=c("G","Adult P"))
RS.plot

ggsave("RS15.pdf",RS.plot,width = 7,height = 5)

#RS 5traits--------
RS_G<-readRDS("RS_G5.Rds")
RS_P<-readRDS("RS_P5.Rds")
RS_P$corr[RS_P$corr>1]<-1
RS_G$corr[RS_G$corr>1]<-1
# show_col(names(sort(table(shade.content[[1]]$raster))))

# shade.content[[3]]$raster[grep("00000000",shade.content[[3]]$raster)]<-colsilh[3]
# cols<-alply(names(sort(table(shade.content[[3]]$raster))),1, color.id)

RS.plot<-
  rbind(RS_G, RS_P) %>% mutate(., Gen=factor(Gen,levels = unique(RS_G$Gen)[c(2:5,1)])) %>%
  subset(., .id!="adults") %>%
  ggplot(aes(Gen,corr))+
  geom_point(aes(size=n, shape=target_type), 
             position=position_jitterdodge(jitter.width = 0.25), alpha=1)+
  # geom_label_repel(aes(label=.id, fill=target_type),position=position_jitterdodge(jitter.width = 0.25))+
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
  scale_shape_manual(values=c(19,21), name="Target Matrix",labels=c("G","Adult P"))+
  scale_fill_manual(values=c("gray","white"), name="Target Matrix",labels=c("G","Adult P"))
RS.plot

ggsave("RS5.pdf",RS.plot,width = 7,height = 5)

#KRZ--------

obsKrz<-readRDS("obsKrz5.Rds")
randKrz<-readRDS("randKrz5.Rds")

krzsubspace.plot<-
  ggplot(melt(obsKrz), aes(as.numeric(variable), value))+
  facet_wrap(.~.id,scales = "free_y",labeller = labeller(.id = c("Monodelphis_B" = "B. Monodelphis (B)","Monodelphis_D" = "A. Monodelphis (D)","Calomys" = "E. Calomys","Cebus"= "D. Sapajus","Didelphis"="C. Didelphis")))+
  # stat_summary(aes(color=.id),fun = median, fun.min = function(x) quantile(x,0.025), fun.max = function(x) quantile(x,0.975))+
  stat_summary(fun = median, fun.min = function(x) quantile(x,0.025), fun.max = function(x) quantile(x,0.975))+
  # scale_x_continuous(limits = c(0,5))+
  geom_line(aes(Var2,value, group=Var1),randKrz, linetype=2)+
  ylab("eigenvalues")+
  xlab("rank")+
  # scale_color_brewer(name="",type = "qual",palette = "Set2")+
  theme(legend.position = "none")
krzsubspace.plot
ggsave("krzsubspace5traits.pdf",krzsubspace.plot,width = 7,height = 5)

obsKrz<-readRDS("obsKrz15.Rds")
randKrz<-readRDS("randKrz15.Rds")

krzsubspace.plot<-
  ggplot(melt(obsKrz), aes(as.numeric(variable), value))+
  facet_wrap(.~.id,scales = "free_y",labeller = labeller(.id = c("Monodelphis_B" = "B. Monodelphis (B)","Monodelphis_D" = "A. Monodelphis (D)","Calomys" = "E. Calomys","Cebus"= "D. Sapajus","Didelphis"="C. Didelphis")))+
  # stat_summary(aes(color=.id),fun = median, fun.min = function(x) quantile(x,0.025), fun.max = function(x) quantile(x,0.975))+
  stat_summary(fun = median, fun.min = function(x) quantile(x,0.025), fun.max = function(x) quantile(x,0.975))+
  # scale_x_continuous(limits = c(0,5))+
  geom_line(aes(Var2,value, group=Var1),randKrz, linetype=2)+
  ylab("eigenvalues")+
  xlab("rank")+
  # scale_color_brewer(name="",type = "qual",palette = "Set2")+
  theme(legend.position = "none")
krzsubspace.plot
ggsave("krzsubspace15traits.pdf",krzsubspace.plot,width = 7,height = 5)
