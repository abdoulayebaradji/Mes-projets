####### Partie 1 : Classification supervisée ####################
library(MASS)
library(ggplot2)
library(GGally)
install.packages("effects")
library(effects)
library(ggeffects)
library(corrplot)
install.packages("forestmodel")
library(forestmodel)

d_cs<-read.csv("/Users/abdoulayebaradji/Downloads/animaux_I.csv",header=TRUE)

summary(d_cs)
# matrice de corrélation 
mat<- cor(d_cs[,-c(1,18)])
corrplot(mat)
# remarque : 
# la variable hair a une forte correlation avec les varibales : eggs milk
# eggs a une forte correlation avec la variable milk
# nous observons de même une fprte cprrélation de la variable milk avec hair et eggs
# backbome et la variable tail aussi ssont corrélées 

d_cs_new<- d_cs[,-c(2,4,10,3,9,11,13)]
mat<- cor(d_cs_new[,-c(1,11)])
corrplot(mat)


summary(d_cs_new)
# après reduction de nos variables qui sont corrélées entre eux pour
# ne pas biaiser l'analyse nous allons maintenant considerer les variables 
# qui sont qualitatives
# 101 observations 18 variables au depart reduit à 11  : dont 10 sont qualitatives 
# la variable classe ici est classe_type
head(d_cs_new)
str(d_cs_new$animal_name)
str(d_cs_new$class_type)


# Recodage des variables factorielles

varName <- names(d_cs_new)[!names(d_cs_new)%in%c("animal_name")]
for (v in varName){
  d_cs_new[,v] <- as.factor(d_cs_new[,v])
}
levels(d_cs_new$class_type) <- c("Mammal","Bird","Reptile","Fish","Amphibian","Bug","Invertebrate")
summary(d_cs_new)
# remarque : nous pouvons voir que les categories de la classe legs présentent des groupes qui sont mal représentés

levels(d_cs_new$legs)[4]<- "8"
#levels(d_cs_new$legs)[5]<-"6"
#levels(d_cs_new$legs)
summary(d_cs_new)

# recodage de la modalité de référence très importantes
# pour l'interprétation des résultats 
# il faut prendre celle qui est bien représentée
summary(d_cs_new)

# Subdivisions échantillons test et échantillons apprentissage

N <- nrow(d_cs_new)
ntrain <- floor(0.75*N)
ntest <- N - ntrain

indices <- sample(1:N,ntrain,replace = FALSE)

d_cs_train <- d_cs_new[indices,]
d_cs_test <- d_cs_new[-indices,]

################## Modèle 2  Forêts aléatoires #########################

# Fonction qui calcule le taux d'erreur
tx_er <- function(pred,vrais){
  mc <- table(pred,vrais)
  1 - sum(diag(mc))/sum(mc)
}

# RF
library(randomForest)
help("randomForest")
mrf <- randomForest(d_cs_train$class_type~.,data=d_cs_train[,-1],method="class",ntree=500,mtry=8)
predrf <- predict(mrf,newdata=d_cs_test,type="class")
te_rf <- tx_er(predrf,d_cs_test$class_type) # 0.15

# Effet du nombre de variables à chaque coupure
effetmtry <- function(m){
  err <- sapply(1:10, FUN = function(i){
    rf <- randomForest(d_cs_train$class_type ~ .,data=d_cs_train,ntree=20,mtry=m)
    predrf <-predict(rf,newdata = d_cs_test)
    return(tx_er(d_cs_test$class_type,predrf))
  })
  return(err)
}

mtryval <- c(1,2,5,8)
err_fn_mtry <- sapply(mtryval,FUN = function(m){effetmtry(m)})
dmBag <- as.data.frame(err_fn_mtry)
names(dmBag) <- paste0("m=",mtryval)
boxplot(dmBag,ylab="Taux d'erreur")

plot(mtryval,apply(err_fn_mtry,2,mean),type="b",xlab="mtry",ylab="Taux d'erreur")

# Effet du nb d'arbres dans la forêt
effetnbTrees <- function(m){
  err <- sapply(1:10, FUN = function(i){
    rf <- randomForest(d_cs_train$class_type ~ .,data=d_cs_train,ntree=m)
    predrf <-predict(rf,newdata = d_cs_test)
    return(tx_er(d_cs_test$class_type,predrf))
  })
  return(err)
}

mval <- c(1,2,5,10,20,50,100,200,500)
err_fn_m <- sapply(mval,FUN = function(m){effetnbTrees(m)})
dm <- as.data.frame(err_fn_m)
names(dm) <- paste0("m=",mval)
boxplot(dm,ylab="Taux d'erreur")

plot(mval,apply(err_fn_m,2,mean),type="b",xlab="Nombre d'arbres",ylab="Taux d'erreur")

mrf$importance
varImpPlot(mrf)
mrf$err.rate # taux d'erreur vectorielle de la prédiction sur les données d'entrée
mrf$predicted

mrf$oob.times
