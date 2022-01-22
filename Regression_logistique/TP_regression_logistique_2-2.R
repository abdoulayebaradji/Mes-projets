# TP régression logistique

# Exo 1 : régression logistique binaire
diab <- read.csv("Documents/Enseignement/Lille/M2 ISN/Analyse de données/TP/tp-ad-m2isn/data/diabetes.csv",header=TRUE)
summary(diab)

library(MASS)
reg.log.bin <- glm(Outcome ~ ., data=diab, family = "binomial")
summary(reg.log.bin)

# odds ratio et leurs intervalles de confiance
exp(coef(reg.log.bin))
exp(confint(reg.log.bin))

tt <- cbind(exp(coef(reg.log.bin)),exp(confint(reg.log.bin)))

# Présentation des résultats
library(broom)
tab.res <- tidy(reg.log.bin, conf.int = TRUE, exponentiate = TRUE)
tab.res

library(gtsummary)
tbl_regression(reg.log.bin, exponentiate = TRUE)

library(forestmodel)
forest_model(reg.log.bin)


# Sélection de modèles
sel.reg <- step(reg.log.bin, k=log(nrow(diab)))
summary(sel.reg)

tbl_regression(sel.reg, exponentiate = TRUE)
forest_model(sel.reg)

# Représentation graphique
library(ggplot2)
ggplot(data=tab.res, aes(x=estimate, y=term, xmin=conf.low, xmax=conf.high)) + geom_point() + geom_errorbar() + geom_vline(xintercept = 1, linetype=2)

library(GGally)
ggcoef(reg.log.bin, exponentiate = TRUE)
ggcoef(reg.log.bin, exponentiate = TRUE, exclude_intercept = TRUE)

library(effects)
library(ggeffects)
ggef <- ggeffect(reg.log.bin, "Glucose")
plot(ggef)
ggef <- ggeffect(reg.log.bin, "Pregnancies")
plot(ggef)

cowplot::plot_grid(plotlist = plot(ggeffect(reg.log.bin)))


# Courbe ROC
pred <- predict(reg.log.bin, type = "response")

sens <- numeric()
spec <- numeric()
for (s in seq(0.001,0.999,0.001)){
  class <- ifelse(pred>s,1,0)
  class <- factor(class,levels=c(0,1))
  mc <- table(class,diab$Outcome)
  sens <- c(sens,mc[1,1]/sum(mc[,1]))
  spec <- c(spec,mc[2,2]/sum(mc[,2]))
}
diffx <- (1-spec[2:999])-(1-spec[1:998])
auc <- sum(sens[2:999]*diffx)

plot(1-spec,sens,type="l",xlim=c(0,1),ylim=c(0,1))
abline(a=0,b=1)

ggplot(data=data.frame(x=1-spec,y=sens),aes(x,y)) + geom_line() + geom_abline(slope=1,intercept = 0, linetype=2) +
  annotate("text",x=0.1,y=0.925,label=paste("AUC=",round(auc,4)))

library(pROC)
roc <- roc(diab$Outcome,pred)
plot(roc)
auc(roc)

# Exo 2 : régression logistique multinomiale
data <- read.table("Documents/Enseignement/Lille/M2 ISN/Analyse de données/TP/tp-ad-m2isn/data/cmc.data",sep=",",header = TRUE)
summary(data)

# Recodage des variables factorielles
varNames <- names(data)[!names(data)%in%c("age","nbchildren")]
for (v in varNames){
  data[,v] <- as.factor(data[,v])
}
levels(data$contraceptive) <- c("Aucune","Court-terme","Long-terme")
summary(data)

# Modèle
library(nnet)
reg.log.mul <- multinom(contraceptive ~ ., data = data)
summary(reg.log.mul)

exp(coef(reg.log.mul))
tab.res <- tidy(reg.log.mul, conf.int = TRUE, exponentiate = TRUE)
tab.res

# Visualisation
ggcoef(reg.log.mul, exponentiate = TRUE) + facet_grid(~y.level)

cowplot::plot_grid(plotlist = plot(ggeffect(reg.log.mul)))


