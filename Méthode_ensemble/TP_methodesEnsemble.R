# TP méthodes d'ensembles

# Importation de la table
d <- read.table("~/Documents/Enseignement/Lille/M2 ISN/AD/TP/Data/newHIV-1_data/1625Data.txt",sep=",",header=T)

head(d)
d$Octamer <- as.character(d$Octamer)
# on sépare les mots de 8 lettres en 8 variables
variables <- t(sapply(1:nrow(d),FUN=function(i){unlist(strsplit(d[i,1],""))}))

d <- cbind(Clived=d$Clived,variables)
summary(d)
d <- as.data.frame(d)
d$Clived <- as.factor(as.numeric(d$Clived)-1)


# Fixe la graine pour obtenir les mêmes résultats à chaque tirage
set.seed(1234)

# Echantillon apprentissage et test
library(caret)
indtrain <- createDataPartition(d$Clived,p=0.8,list=F)
dtrain <- d[indtrain,]
dtest <- d[-indtrain,]

# Proportion d'octamères clivés dans les deux bases
prop.table(table(dtrain$Clived))
prop.table(table(dtest$Clived))

# Fonction qui calcule le taux d'erreur
tx_er <- function(pred,vrais){
  mc <- table(pred,vrais)
  1 - sum(diag(mc))/sum(mc)
}

### Classifieur constant
nbClived <- table(dtrain$Clived)[2]
mcst <- mean(dtest$Clived!=0)

### CART
# Arbre max
library(rpart)
mcartmax <- rpart(as.factor(Clived)~.,data=dtrain,cp=0,minbucket=1,maxdepth=30)
predcartmax <- predict(mcartmax,newdata=dtest,type="class")
te_cartmax <- tx_er(predcartmax,dtest$Clived) 

# Elagage
mcart <- prune(mcartmax,cp=0.021)
predcart <- predict(mcart,newdata=dtest,type="class")
te_cart <- tx_er(predcart,dtest$Clived)

# Arbre à 1 noeud
mdecstump <- rpart(as.factor(Clived)~.,data=dtrain,cp=0,maxdepth=1)
preddecstump <- predict(mdecstump,newdata=dtest,type="class")
te_decst <- tx_er(preddecstump,dtest$Clived) 

### Bagging
library(adabag)
# Par défaut
mbag <- bagging(Clived~.,data=dtrain,mfinal=20)
predbag20 <- predict(mbag,newdata=dtest,type="class")
te_bag20 <- tx_er(predbag20$class,dtest$Clived) 

# Bagging d'arbres à 1 noeud
mbagstump <- bagging(Clived~.,data=dtrain,mfinal=20,control=rpart.control(cp=0,maxdepth=1,minbucket=1))
predbagstump <- predict(mbagstump,newdata=dtest,type="class")
te_bagstump <- tx_er(predbagstump$class,dtest$Clived) 

# Bagging d'arbres profonds
mbagdeep <- bagging(Clived~.,data=dtrain,mfinal=20,control=rpart.control(cp=0,maxdepth=30,minbucket=1))
predbagdeep <- predict(mbagdeep,newdata=dtest,type="class")
te_bagdeep <- tx_er(predbagdeep$class,dtest$Clived)


# Visualisation de l'effet du nombre d'arbres
effetnbTreesBag <- function(m){
  err <- sapply(1:10, FUN = function(i){
    bag <- bagging(Clived ~ .,data=dtrain,mfinal=m,control=rpart.control(cp=0,maxdepth=30,minbucket=1))
    predbag <-predict(bag,newdata = dtest)
    return(tx_er(dtest$Clived,predbag$class))
  })
  return(err)
}

mval <- c(1,2,5,10,20,50)
err_fn_m_bag <- sapply(mval,FUN = function(m){effetnbTreesBag(m)})
dmBag <- as.data.frame(err_fn_m_bag)
names(dmBag) <- paste0("m=",mval)
boxplot(dmBag,ylab="Taux d'erreur")

plot(mval,apply(err_fn_m_bag,2,mean),type="b",xlab="Nombre d'arbres baggés", ylab="Taux d'erreur")



# RF
library(randomForest)
mrf <- randomForest(Clived~.,data=dtrain,method="class",ntree=20,mtry=1)
predrf <- predict(mrf,newdata=dtest,type="class")
te_rf <- tx_er(predrf,dtest$Clived) 

# Effet du nombre de variables à chaque coupure
effetmtry <- function(m){
  err <- sapply(1:10, FUN = function(i){
    rf <- randomForest(Clived ~ .,data=dtrain,ntree=20,mtry=m)
    predrf <-predict(rf,newdata = dtest)
    return(tx_er(dtest$Clived,predrf))
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
    rf <- randomForest(Clived ~ .,data=dtrain,ntree=m)
    predrf <-predict(rf,newdata = dtest)
    return(tx_er(dtest$Clived,predrf))
  })
  return(err)
}

mval <- c(1,2,5,10,20,50,100,200,500)
err_fn_m <- sapply(mval,FUN = function(m){effetnbTrees(m)})
dm <- as.data.frame(err_fn_m)
names(dm) <- paste0("m=",mval)
boxplot(dm,ylab="Taux d'erreur")

plot(mval,apply(err_fn_m,2,mean),type="b",xlab="Nombre d'arbres",ylab="Taux d'erreur")


# Boosting
mbooststump <- boosting(Clived~.,data=dtrain,mfinal=20,boos=F,control=rpart.control(cp=0,maxdepth=1,minbucket=1))
predbooststump <- predict(mbooststump,newdata=dtest,type="class")
te_booststump <- tx_er(predbooststump$class,dtest$Clived) 

mboostdeep <- boosting(Clived~.,data=dtrain,mfinal=20,boos=F,control=rpart.control(cp=0,maxdepth=30,minbucket=1))
predboostdeep <- predict(mboostdeep,newdata=dtest,type="class")
te_boostdeep <- tx_er(predboostdeep$class,dtest$Clived) 


# XGBoost
library(xgboost)
library(ade4)

# Transforme les données en dummy variables
class <- as.numeric(d$Clived)
class <- class-1
ddummy <- cbind(class,acm.disjonctif(d[,-1]))

ddumtrain <- ddummy[indtrain,]
ddumtest <- ddummy[-indtrain,]


dtrainXG <- xgb.DMatrix(as.matrix(ddumtrain[,-1]),label=as.matrix(ddumtrain[,1]))
dtestXG <- xgb.DMatrix(as.matrix(ddumtest[,-1]),label=as.matrix(ddumtest[,1]))
watchlist <- list(train=dtrainXG,test=dtestXG)

mxgb <- xgb.train(params = list(max_depth=30, eta = 0.75, objective="binary:logistic"),
                  data = dtrainXG,
                  nrounds = 20,
                  watchlist = watchlist)
