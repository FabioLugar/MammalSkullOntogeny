setwd("~/Onto1sp_bayesian/")
library(plyr)
library(dplyr)
library(ggplot2)
library(RRPP)
library(cowplot)
library(ggforce)
models<-readRDS(file = "models.Rds")

data<-models %>% llply(., function(x){
  ldply(x, function(y)  {
    if(is.null(y$model$`as.factor(ESPECIE)`)) sp=NA else sp= as.factor(y$model$`as.factor(ESPECIE)`)
    if(is.null(y$model$`as.factor(SEXO)`)) sex=NA else sex= as.factor(y$model$`as.factor(SEXO)`)
    data.frame(sp, sex, y$model$temp.Y)
  },.inform = T)
},.inform = T)

cebusdata<-readRDS("CEBUSData.RDS")
data$Cebus$sp<-cebusdata$ESPECIE
data$Cebus$sp[data$Cebus$sp=="sp"]<-"sp."
x<-data$Cebus
x<-subset(x, .id!="adult")
x$.id<-factor(x$.id, unique(x$.id))
pca<-prcomp(x[,-c(1:3)],scale. = F)
pca_cebus.plot<-ggplot(data.frame(x[,1:3],pca$x), aes(PC1,PC2))+
  geom_point(aes(color=sp))+
  coord_fixed()+
  scale_color_brewer(name="Species",type = "qual",palette=2)+
  theme_classic()+
  scale_y_continuous(expand = expansion(add = 6))

rrpp.df<-rrpp.data.frame(vars=as.matrix(x[,-c(1:3)]),size=rowSums(log(x[,-c(1:3)])),sp=x$sp)
lmrrpp3<-lm.rrpp(vars~size*sp,data = rrpp.df,Parallel = T,SS.type="II")
anova(lmrrpp3)


ggsave("pca_cebus.pdf",pca_cebus.plot,width = 6,height = 3)

write.csv(anova(lmrrpp3) $ table,file = "npmanova_cebus.csv")


#############

load("raw.data.RData")
info_virginiana<-read.csv("dados_de.csv")

domestica<-raw.data.no.outlier.only.good.sample$Didelphimorphia %>% subset(., ESPECIE=="domestica")
virginiana<-raw.data.no.outlier.only.good.sample$Didelphimorphia %>% subset(., ESPECIE=="virginiana")
names(info_virginiana)[1:5]<-names(virginiana)[1:5]

virginiana_info<-left_join(virginiana,info_virginiana)
dim(virginiana_info)
pca<-select(virginiana_info,ISPM:BAOPI) %>% prcomp

table(virginiana_info$Estado)

cbind(virginiana_info, pca$x) %>%
  filter(., !is.na(Estado)) %>% 
  ggplot(., aes(PC1,PC2))+
  theme_classic()+
  geom_point(aes(shape=as.factor(AGE), color=Estado), show.legend = T)+
  scale_shape_manual(values=1:8)+
  coord_fixed()+
  scale_y_continuous(expand = expansion(add = 20))+
  stat_ellipse(aes(color=Estado))
  