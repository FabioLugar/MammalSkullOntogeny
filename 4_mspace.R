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
registerDoParallel(cores=50)
matricesDist<-readRDS("matricesDist.Rds")
info<-readRDS("ids.Rds")
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
saveRDS(prcoa,"prcoa.Rds")


multipcoa<-foreach(i=levels(info$gen)) %do% {
  I<-info$gen==i & info$age!="adults"
  mats<-matricesDist[I] %>%
    llply(., function(x) {
      x<-x$P
      x/psych::tr(x)
    })
    
  matricesA<-abind(mats, along = 3)
  meanP<-apply(matricesA,1:2,weighted.mean,w=info$ns[I])
  eigenP<-eigen(meanP)
  
  vargrad<-ExtendMatrix(meanP)$var.grad
  dims<-vargrad$PC[which(vargrad$Gradient.Variance<1.e-07)[1]]
  dims<-20
  matricesLs<-matricesDist[I] %>%
    llply(., function(x) {
      for(i in 1:100) {
        X<-t(eigenP$vectors[,1:dims]) %*% x$Ps[i,,] %*% eigenP$vectors[,1:dims]
        x$Ps[i,1:dims,1:dims]<-X/psych::tr(X)
      }
      alply(x$Ps[1:100,1:dims,1:dims], 1, extract)
    })
  matricesLs<-unlist(matricesLs, recursive = FALSE)
  eigen.phen<-MatrixDistance(matricesLs, distance = c("RiemannDist"),parallel = TRUE)
  prcoa <- pr.coord(as.matrix(as.dist(eigen.phen)))
}
names(multipcoa)<-levels(info$gen)
saveRDS(multipcoa,"multiprcoa.Rds")
