models<-readRDS(file = "models.Rds")

temp.tdd<-subset(models$Monodelphis_B$adults$model,`as.factor(AGE)`%in% c(300,400))
names(temp.tdd)[2:3]<-c("SEXO", "AGE")
models$Monodelphis_B$adults<-lm(formula = temp.Y ~ as.factor(SEXO) + as.factor(AGE), data = temp.tdd)

dim(models$Monodelphis_B$adults$model)
saveRDS(models,"models.Rds")
