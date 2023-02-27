
repKC<-c(llply(models$Monodelphis_D, function(x){
  BootstrapRep(x$residuals, KrzCor, iterations = 100, parallel = T)
}),
llply(models$Monodelphis_B, function(x){
  BootstrapRep(x$residuals, KrzCor, iterations = 100, parallel = T)
}),
llply(models$Didelphis,  function(x){
  BootstrapRep(x$residuals, KrzCor, iterations = 100, parallel = T)
}),
llply(models$Cebus,  function(x){
  BootstrapRep(x$residuals, KrzCor, iterations = 100, parallel = T)
}),
llply(models$Calomys, function(x){
  BootstrapRep(x$residuals, KrzCor, iterations = 100, parallel = T)
}))
repKC<-unlist(repKC)
