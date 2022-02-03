install.packages("ggplot2")
library(FactoMineR)
library(factoextra)

#############Partie 1 :Représentation graphique de nos simulations #################### 
donne0<-read.table("/Users/abdoulayebaradji/Desktop/cranck_nicolson.txt",sep=" ",header=FALSE)
donne0_wuhan<-read.table("/Users/abdoulayebaradji/Desktop/cranck_nicolson_wuhan.txt",sep=" ",header=FALSE)
donne1<-read.table("/Users/abdoulayebaradji/Desktop/heun.txt",sep=" ",header=FALSE)
donne1_wuhan<-read.table("/Users/abdoulayebaradji/Desktop/heun_wuhan.txt",sep=" ",header=FALSE)
donne2<-read.table("/Users/abdoulayebaradji/Desktop/RK4.txt",sep=" ",header=FALSE)
donne2_wuhan<-read.table("/Users/abdoulayebaradji/Desktop/RK4_wuhan.txt",sep=" ",header=FALSE)
donne0_modified<-read.table("/Users/abdoulayebaradji/Desktop/cranck_modified.txt",sep=" ", header=FALSE)
donne1_modified<-read.table("/Users/abdoulayebaradji/Desktop/heun_modified.txt",sep=" ",header=FALSE)
donne2_modified<-read.table("/Users/abdoulayebaradji/Desktop/RK4_modified.txt",sep=" ",header=FALSE)
donne<-read.table('/Users/abdoulayebaradji/Desktop/erreur_cranck.txt',sep=" ",header=FALSE)
donne_modified<-read.table('/Users/abdoulayebaradji/Desktop/erreur_cranck_modified.txt',sep=" ",header=FALSE)
donne3<-read.table('/Users/abdoulayebaradji/Desktop/erreur_heun.txt',sep=" ",header=FALSE)
donne3_modified<-read.table('/Users/abdoulayebaradji/Desktop/erreur_heun_modified.txt',sep=" ",header=FALSE)
donne4<-read.table('/Users/abdoulayebaradji/Desktop/erreur_rk4.txt',sep=" ",header=FALSE)
donne4_modified<-read.table('/Users/abdoulayebaradji/Desktop/erreur_rk4_modified.txt',sep=" ",header=FALSE)
donne_wuhan<-read.table('/Users/abdoulayebaradji/Desktop/erreur_cranck_wuhan.txt',sep=" ",header=FALSE)
donne3_wuhan<-read.table('/Users/abdoulayebaradji/Desktop/erreur_heun_wuhan.txt',sep=" ",header=FALSE)
donne4_wuhan<-read.table('/Users/abdoulayebaradji/Desktop/erreur_rk4_wuhan.txt',sep=" ",header=FALSE)

############################## schéma SEIR population N=1########################## 
plot(seq(0,30,length=1025),log(donne0[,1]),type="l",main="evolution de la population aucours du temps schéma cranck")
lines(seq(0,30,length=1025),log(donne0[,2]),type="l",col="red")
lines(seq(0,30,length=1025),log(donne0[,3]),type="l",col="blue")
lines(seq(0,30,length=1025),log(donne0[,4]),type="l",col="green")
legend("topright",legend=c("Sane","Exposed","Infected", "Recovered"),
       col=c("black", "red","blue","green"), lty=1:2, cex=0.8)
plot(seq(0,30,length=1025),log10(apply(donne0,1,sum)),main="population finale totale",type="l")

plot(seq(0,90,length=925),log10(donne0_wuhan[,1]),type="l",main="évolution population de wuhan au cours du temps schéma cranck")
lines(seq(0,90,length=925),log10(donne0_wuhan[,2]),type="l",col="red")
lines(seq(0,90,length=925),log10(donne0_wuhan[,3]),type="l",col="blue")
lines(seq(0,90,length=925),log10(donne0_wuhan[,4]),type="l",col="green")
legend("topright",legend=c("Sane","Exposed","Infected", "Recovered"),
       col=c("black", "red","blue","green"), lty=1:2, cex=0.8)
plot(seq(0,90,length=925),(apply(donne0_wuhan,1,sum)),main="population finale totale wuhan schéma cranck",type="l")

plot(seq(0,30,length=1001),log(donne1[,1]),type="l",main="évolution population aucours du temps schéma heun")
lines(seq(0,30,length=1001),log(donne1[,2]),type="l",col="red")
lines(seq(0,30,length=1001),log(donne1[,3]),type="l",col="blue")
lines(seq(0,30,length=1001),log(donne1[,4]),type="l",col="green")
legend("topright",legend=c("Sane","Exposed","Infected", "Recovered"),
       col=c("black", "red","blue","green"), lty=1:2, cex=0.8)
plot(seq(0,30,length=1001),log(apply(donne1,1,sum)),main="population finale totale",type="l")

################# schéma SEIR population de wuhan #######################################
plot(seq(0,90,length=901),log10(donne1_wuhan[,1]),type="l",main="évolution population wuhan aucours du temps chéma heun")
lines(seq(0,90,length=901),log10(donne1_wuhan[,2]),type="l",col="red")
lines(seq(0,90,length=901),log10(donne1_wuhan[,3]),type="l",col="blue")
lines(seq(0,90,length=901),log10(donne1_wuhan[,4]),type="l",col="green")
legend("topright",legend=c("Sane","Exposed","Infected", "Recovered"),
       col=c("black", "red","blue","green"), lty=1:2, cex=0.8)
plot(seq(0,90,length=901),(apply(donne1_wuhan,1,sum)),main="population finale totale wuhan schema heun",type="l")


plot(seq(0,30,length=1001),log(donne2[,1]),type="l",main="évolution population aucours du temps schema RK4")
lines(seq(0,30,length=1001),log(donne2[,2]),type="l",col="red")
lines(seq(0,30,length=1001),log(donne2[,3]),type="l",col="blue")
lines(seq(0,30,length=1001),log(donne2[,4]),type="l",col="green")
legend("topright",legend=c("Sane","Exposed","Infected", "Recovered"),
       col=c("black", "red","blue","green"), lty=1:2, cex=0.8)

plot(seq(0,90,length=901),log10(donne2_wuhan[,1]),type="l",main="évolution population wuhan aucours du temps schéma RK4")
lines(seq(0,90,length=901),log10(donne2_wuhan[,2]),type="l",col="red")
lines(seq(0,90,length=901),log10(donne2_wuhan[,3]),type="l",col="blue")
lines(seq(0,90,length=901),log10(donne2_wuhan[,4]),type="l",col="green")
legend("topright",legend=c("Sane","Exposed","Infected", "Recovered"),
       col=c("black", "red","blue","green"), lty=1:2, cex=0.8)
plot(seq(0,90,length=901),log(apply(donne2_wuhan,1,sum)),main="population finale totale wuhan schéma RK4",type="l")

########################## schéma SEIR Modified graphe ##############################
plot(seq(0,60,length=1025),log10(donne0_modified[,1]),type="l",main="Susceptible schéma SEIR Modified")
plot(seq(0,60,length=1025),log10(donne0_modified[,2]),type="l",col="red",main="exposed shcéma SEIR Modified")
plot(seq(0,60,length=1025),log10(donne0_modified[,3]),type="l",col="blue",main="infected schéma SEIR Modified")
plot(seq(0,60,length=1025),log10(donne0_modified[,4]),type="l",col="green",main="asymptomatic schéma SEIR Modified")
plot(seq(0,60,length=1025),log10(donne0_modified[,5]),type="l",col="yellow",main="quarantined susceptible schéma SEIR Modified")
plot(seq(0,60,length=1025),log10(donne0_modified[,6]),type="l",col="orange",main="isolated exposed schéma SEIR Modifiied")
plot(seq(0,60,length=1025),log10(donne0_modified[,7]),type="l",col="brown",main="hospitalized schéma SEIR Modified")
plot(seq(0,60,length=1025),log10(donne0_modified[,8]),type="l",col="purple",main="recovered schéma SEIR Modified")
plot(seq(0,60,length=1001),(apply(donne1_modified,1,sum)),main="variation population finale totale SEIR modified",type="l")


plot(seq(0,60,length=1001),log10(donne1_modified[,8]),type="l",col="purple")
plot(seq(0,60,length=1001),log10(donne1_modified[,1]),type="l",main="évolution population aucours du temps schéma heun modified ")
plot(seq(0,60,length=1001),log10(donne1_modified[,2]),type="l",col="red")
plot(seq(0,60,length=1001),log10(donne1_modified[,3]),type="l",col="blue")
plot(seq(0,60,length=1001),log10(donne1_modified[,4]),type="l",col="green")
plot(seq(0,60,length=1001),log10(donne1_modified[,5]),type="l",col="yellow")
plot(seq(0,60,length=1001),log10(donne1_modified[,6]),type="l",col="orange")
plot(seq(0,60,length=1001),log10(donne1_modified[,7]),type="l",col="brown")


plot(seq(0,60,length=1001),log10(donne2_modified[,1]),type="l",main="évolution population aucours du temps schéma rk4 modified")
plot(seq(0,60,length=1001),log10(donne2_modified[,2]),type="l",col="red")
plot(seq(0,60,length=1001),log10(donne2_modified[,3]),type="l",col="blue")
plot(seq(0,60,length=1001),log10(donne2_modified[,4]),type="l",col="green")
plot(seq(0,60,length=1001),log10(donne2_modified[,5]),type="l",col="yellow")
plot(seq(0,60,length=1001),log10(donne2_modified[,6]),type="l",col="orange")
plot(seq(0,60,length=1001),log10(donne2_modified[,7]),type="l",col="brown")
plot(seq(0,60,length=1001),log10(donne2_modified[,8]),type="l",col="purple")

########################## graphique des erreurs :################################## 
h<-numeric(10)
h[1]<-1
for (i in 1:9){
  h[i+1]<-1/2**i
}

plot(log10(h),log10(donne[,1]),type="o",main="erreur modèle SEIR population N=1 schema cranck")
abline(0,2,col="red")
abline(0,4,col="blue")

plot(log10(h),log10(donne3[,1]),type="o",main="erreur modèle SEIR population N=1 schéma heun")
abline(0,2,col="red")
abline(0,4,col="blue")

plot(log10(h),log10(donne4[,1]),type="o",main="erreur modèle SEIR population N=1 schéma rk4")
plot(log10(h),log10(donne_wuhan[,1]),type="o",main="erreur modèle SEIR population wuhan schéma cranck")
plot(log10(h),log10(donne3_wuhan[,1]),type="o",col="red",main="erreur modèle SEIR population wuhan schéma heun")
abline(0,2,col="red")
abline(0,4,col="blue")

plot(log10(h),log10(donne4_wuhan[,1]),type="o",col="blue",main="erreur modèle SEIR population wuhan schéma rk4")
abline(0,2,col="red")
abline(0,4,col="blue")
plot(log10(h),log10(donne_modified[,1]),type="o",main="erreur pour le modèle SEIR modifed schéma cranck")
plot(log10(h),log10(donne3_modified[,1]),type="o",col="red",main="erreur modèle SEIR modified schéma heun")
plot(log10(h),log10(donne4_modified[,1]),type="o",col="blue",main="erreur modèle SEIR modified schéma rk4")


#################### Partie 2 : analyse de donnée ###########################

data0<- read.csv("/Users/abdoulayebaradji/Desktop/covid_confirmed_deaths_since_5th.csv",header=TRUE)
data1<- read.csv("/Users/abdoulayebaradji/Desktop/covid_confirmed_cases_since_100th.csv",header=TRUE)
data0_death_franc<-data0[7068:7214,]
data0_death_ital<-data0[10090:10236,]
data0_death_spai<-data0[18267:18412,]

# l'évolution de la trajectoire de chaque pays à partir 5 décès : 

data0_death_france<-data0_death_franc[which(data0_death_franc$X.deaths. > 4),]
head(data0_death_france)
summary(data0_death_france)

data0_death_italy<-data0_death_ital[which(data0_death_ital$X.deaths. > 4),]
head(data0_death_italy)
summary(data0_death_italy)

data0_death_spain<-data0_death_spai[which(data0_death_spai$X.deaths. > 4),]
head(data0_death_spain)
summary(data0_death_spain)

plot(data0_death_france$Days.since.the.5th.total.confirmed.death,log(data0_death_france$X.deaths.),type="l",col="black",xlab="Days since the 5th confirmed death",main="vitesse à laquelle le nombre de décès confirmés a augmenté ")
lines(data0_death_italy$Days.since.the.5th.total.confirmed.death,log(data0_death_italy$X.deaths.),type="l",col="red")
lines(data0_death_spain$Days.since.the.5th.total.confirmed.death,log(data0_death_spain$X.deaths.),type="l",col="blue")
legend("topleft",legend=c("death France", "death Italy","death Spain"),
       col=c("black", "red","blue"), lty=1:2, cex=0.8)
# en prennant le logarithme : 

plot(data0_death_france$Days.since.the.5th.total.confirmed.death,log(data0_death_france$X.deaths.),xlab="Days since the 100th confirmed death",ylab="LOG",type="l",col="black",main="logarithme vitesse à laquelle le nombre décès confirmés a augmenté")
lines(data0_death_italy$Days.since.the.5th.total.confirmed.death,log(data0_death_italy$X.deaths.),type="l",col="red")
lines(data0_death_spain$Days.since.the.5th.total.confirmed.death,log(data0_death_spain$X.deaths.),type="l",col="blue")
legend("topleft",legend=c("death France", "death Italy","death Spain"),
       col=c("black", "red","blue"), lty=1:2, cex=0.8)

# l'evolution de vitesse à laquelle le nombre de cas a augmenté dans chque pays après 100 cas : 

data1_case_franc<-data1[7088:7234,]
data1_case_ital<-data1[10110:10256,]
data1_case_spai<-data1[18273:18418,]

data1_case_france<-data1_case_franc[which(data1_case_franc$X.cases. > 99),]
head(data1_case_france)
summary(data1_case_france)

data1_case_italy<-data1_case_ital[which(data1_case_ital$X.cases. > 99),]
head(data1_case_italy)
summary(data1_case_italy)

data1_case_spain<-data1_case_spai[which(data1_case_spai$X.cases. > 99),]
head(data1_case_spain)
summary(data1_case_spain)

plot(data1_case_france$Days.since.the.100th.confirmed.case..days.,(data1_case_france$X.cases.),type="l",col="black",xlab="Days since the 100th confirmed cases",main="vitesse à laquelle le nombre de cas confirmés a augmenté ")
lines(data1_case_italy$Days.since.the.100th.confirmed.case..days.,(data1_case_italy$X.cases.),type="l",col="red")
lines(data1_case_spain$Days.since.the.100th.confirmed.case..days.,(data1_case_spain$X.cases.),type="l",col="blue")
legend("topleft",legend=c("infected France", "infected Italy","infected Spain"),
       col=c("black", "red","blue"), lty=1:2, cex=0.8)
# en prennant le logarithme : 

plot(data1_case_france$Days.since.the.100th.confirmed.case..days.,log10(data1_case_france$X.cases.),xlab="Days since the 100th confirmed cases",ylab="LOG",type="l",col="black",main="logarithme vitesse à laquelle le nombre cas confirmés a augmenté")
lines(data1_case_italy$Days.since.the.100th.confirmed.case..days.,log10(data1_case_italy$X.cases.),type="l",col="red")
lines(data1_case_spain$Days.since.the.100th.confirmed.case..days.,log10(data1_case_spain$X.cases.),type="l",col="blue")
legend("topleft",legend=c("infected France", "infected Italy","infected Spain"),
       col=c("black", "red","blue"), lty=1:2, cex=0.8)

# NOS SIMULATIONS PAR LE MODÈLE SEIR ET COMPARAISON AVEC LES DONNÉES RÉELLES : 

data2<-read.table("/Users/abdoulayebaradji/Desktop/RK4_france.txt",sep=" ",header=FALSE)
colnames(data2)<-c("S","E","I","R")
data3<-read.table("/Users/abdoulayebaradji/Desktop/RK4_italy.txt",sep=" ",header=FALSE)
colnames(data3)<-c("S","E","I","R")
data4<-read.table("/Users/abdoulayebaradji/Desktop/RK4_spain.txt",sep=" ",header=FALSE)
colnames(data4)<-c("S","E","I","R")
data2_isoleted<-read.table("/Users/abdoulayebaradji/Desktop/RK4_france_isolement.txt",sep=" ",header=FALSE)
colnames(data2_isoleted)<-c("S","E","I","R")
data3_isoleted<-read.table("/Users/abdoulayebaradji/Desktop/RK4_italy_isolement.txt",sep=" ",header=FALSE)
colnames(data3_isoleted)<-c("S","E","I","R")
data4_isoleted<-read.table("/Users/abdoulayebaradji/Desktop/RK4_spain_isolement.txt",sep=" ",header=FALSE)
colnames(data4_isoleted)<-c("S","E","I","R")

data2_isoleted<-data2_isoleted[130:901,]
data3_isoleted<-data3_isoleted[130:901,]
data4_isoleted<-data4_isoleted[120:901,]

Ndeath_f<- 67e+6 - apply(data2,1,sum)
Ndeath_i<- 60e+6 - apply(data3,1,sum)
Ndeath_s<- 47e+6 - apply(data4,1,sum)

Ndeath_f_isoleted<- 67e+6 - apply(data2_isoleted,1,sum)
Ndeath_i_isoleted<- 60e+6 - apply(data3_isoleted,1,sum)
Ndeath_s_isoleted<- 47e+6 - apply(data4_isoleted,1,sum)


plot(seq(0,90,length=901),log10(data2$I),type="l",main="France")
lines(seq(0,90,length=86),log10(data1_case_france$X.cases.),type="l",col="red")
lines(seq(0,90,length=901),log10(Ndeath_f),type="l",col="blue")
lines(seq(0,90,length=81),log10(data0_death_france$X.deaths.),type="l",col="green")
lines(seq(0,90,length=772),log10(data2_isoleted$I),type="l",col="brown")
lines(seq(0,90,length=772),log10(Ndeath_f_isoleted),type="l",col="orange")
legend("topleft",legend=c("infected experiment","infected data real","death experiment", "death data real","infected experiement isolation","death experiment isolation"),
       col=c("black", "red","blue","green","brown","orange"), lty=1:2, cex=0.8)

plot(seq(0,90,length=901),log10(data3$I),type="l",main="Italy")
lines(seq(0,90,length=92),log10(data1_case_italy$X.cases.),type="l",col="red")
lines(seq(0,90,length=901),log10(Ndeath_i),type="l",col="blue")
lines(seq(0,90,length=91),log10(data0_death_italy$X.deaths.),type="l",col="green")
lines(seq(0,90,length=772),log10(data3_isoleted$I),type="l",col="brown")
lines(seq(0,90,length=772),log10(Ndeath_i_isoleted),type="l",col="orange")
legend("topleft",legend=c("infected experiment","infected data real","death experiment", "death data real","infected experiment isolation","death experiment isolation"),
       col=c("black", "red","blue","green"), lty=1:2, cex=0.8)

plot(seq(0,90,length=901),log10(data4$I),type="l",main="Spain")
lines(seq(0,90,length=84),log10(data1_case_spain$X.cases.),type="l",col="red")
lines(seq(0,90,length=901),log10(Ndeath_s),type="l",col="blue")
lines(seq(0,90,length=79),log10(data0_death_spain$X.deaths.),type="l",col="green")
lines(seq(0,90,length=782),log10(data4_isoleted$I),type="l",col="brown")
lines(seq(0,90,length=782),log10(Ndeath_s_isoleted),type="l",col="orange")
legend("topleft",legend=c("infected experiment", "infected data real","death experiment","death data real","infected experiment isolation","death experiment isolation"),
       col=c("black", "red","blue","green"), lty=1:2, cex=0.8)

