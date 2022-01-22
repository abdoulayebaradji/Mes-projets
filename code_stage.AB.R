############### importation des packages ##################################
library(ggplot2)
library(readxl)
library(FactoMineR) # MCA
library(factoextra) # visualisation des graphes
library(missMDA) # données manquantes (imputation)
library(tidyverse) # récodage des catégories
library(naniar)
library(corrplot)
library(GGally) 
install.packages("here")
library(here)
################ Importation base de donnée Projet DIFUS_E originale ####################@################
enquete<- read_excel("data/Enquete.xlsx") %>% 
  as.data.frame()

# subdivision de la base de données en fonction des analyses à effectuer
enquete1<- enquete[,-c(2,5,6,7,8,9)]
# on enlève toutes les questions ouvertes
enquete2<- enquete1[,-c(5,7,9,11,13,15,17,19,21,23,25,27,29,31,42,48)]
# on enlève la première ligne de la table
enquete3<- enquete2[-1,]

# renommer les variables
# Q1.? concerne les  questions sur les avis des pratiques collectifs
# Q2.? concerne les questions sur les avis des pratiques individuels
# Q3.? les questions concernant les avis sur la position du CIRAD
# Q4.? concerne les questions socio-démographiques
names(enquete3)<- c("Respondent.ID", "Start.Date", "End.Date", "Q1.1", 
                    "Q1.2", 
                    "Q2.1", 
                    "Q2.2", 
                    "Q2.3", 
                    "Q2.4", 
                    "Q2.5", 
                    "Q2.6", 
                    "Q2.7", 
                    "Q2.8", 
                    "Q2.9", 
                    "Q3.1", 
                    "Q3.2", 
                    "Q3.3", 
                    "Q3.4", 
                    "Q3.5", "Q3.6", "Q3.7", "Q3.8", "Q3.9", "Q3.10", "Q3.11", 
                    "Q3.12", "Q3.13", "Q4.1", "Q4.2", 
                    "Q4.3", "Q4.4", "Q4.5"
)

# On enlève les variables qui ne décrivent pas une opinion etc ...
enquete4<- enquete3[,-c(1,2,3,28,29,30,31,32)]

# 28 : 32 et 1 : 3 seront considérées comme suplémentaires





############ quelques statistiques descriptives sur les variables concernées #############################
#barplot(prop.table(table(enquete$NA..4)))
help("barplot")

# L'urgence climatique va exiger des changements profonds dans la pratique de nos métiers de recherche


barplot(table(enquete3$Q1.1))
barplot(prop.table(table(enquete3.new$Q1.1.new))) 
# L'organisation, l'évaluation et le financement de la recherche poussent à des déplacements trop nombreux et inutiles
barplot(prop.table(table(enquete3$Q1.2)))
barplot(prop.table(table(enquete3.new$Q1.2.new)))

# questions : Avis sur les pratiques individuels
barplot(prop.table(table(enquete3$Q2.1)))
barplot(prop.table(table(enquete3.new$Q2.1.new)))


barplot(prop.table(table(enquete3$Q2.2)))
barplot(prop.table(table(enquete3.new$Q2.2.new)))

barplot(prop.table(table(enquete3$Q2.3)))
barplot(prop.table(table(enquete3$Q2.3.new)))

barplot(prop.table(table(enquete3$Q2.4)))
barplot(prop.table(table(enquete3$Q2.4.new)))

barplot(prop.table(table(enquete3$Q2.5)))
#barplot(prop.table(table(enquete3$Q2.5.new)))

barplot(prop.table(table(enquete3$Q2.6)))
barplot(prop.table(table(enquete3$Q2.6.new)))

barplot(prop.table(table(enquete3$Q2.7)))
barplot(prop.table(table(enquete3$Q2.7.new)))

barplot(prop.table(table(enquete3$Q2.8)))
barplot(prop.table(table(enquete3$Q2.8.new)))

barplot(prop.table(table(enquete3$Q2.9)))
barplot(prop.table(table(enquete3$Q2.9.new)))

# Questions : Avis sur la positions du CIRAD 
barplot(prop.table(table(enquete3$Q3.1)))
barplot(prop.table(table(enquete3.new$Q3.1.new)))

barplot(prop.table(table(enquete3$Q3.2)))
barplot(prop.table(table(enquete3$Q3.2.new)))

barplot(prop.table(table(enquete3$Q3.3)))
barplot(prop.table(table(enquete3$Q3.3.new)))

barplot(prop.table(table(enquete3$Q3.4)))
barplot(prop.table(table(enquete3.new$Q3.4.new)))

barplot(prop.table(table(enquete3$Q3.5)))
barplot(prop.table(table(enquete3$Q3.5.new)))

barplot(prop.table(table(enquete3$Q3.6)))
barplot(prop.table(table(enquete3$Q3.6.new)))

barplot(prop.table(table(enquete3$Q3.7)))
barplot(prop.table(table(enquete3$Q3.7.new)))

barplot(prop.table(table(enquete3$Q3.8)))
barplot(prop.table(table(enquete3$Q3.8.new)))

barplot(prop.table(table(enquete3$Q3.9)))
barplot(prop.table(table(enquete3$Q3.9.new)))

barplot(prop.table(table(enquete3$Q3.10)))

barplot(prop.table(table(enquete3$Q3.11)))

barplot(prop.table(table(enquete3$Q3.12)))

barplot(prop.table(table(enquete3$Q3.13)))

# Questions : socio-démographiques
barplot(prop.table(table(enquete3.new$Q4.3.new)))

barplot(prop.table(table(enquete3$Q4.2)))

barplot(prop.table(table(enquete3$Q4.3)))
barplot(prop.table(table(enquete3.new$Q4.3.new)))

barplot(prop.table(table(enquete3$Q4.4)))

barplot(prop.table(table(enquete3$Q4.5)))


d_224<- read_excel(here::here("data/tab_224.xlsx"), col_names = TRUE, skip = 1) %>% 
  as.data.frame()


######### ACM sur la base de données enquete 4 en elevant les données manquantes #####################

help(which)
# traitement de données manquantes
is.na(enquete4) # nombre de données manquantes
which(is.na(enquete4),arr.ind = TRUE) # couple ligne/colonne des individus presentant une valeur manquante
# on enlève de la base les individus ayant au moins un NA
indligneNA <- which(is.na(enquete4),arr.ind = TRUE)[,1] # indice des lignes dans la colonne 1
length(unique(indligneNA)) # 69 individus ont au moins une donnée manquante

table(indligneNA) # nombre d'occurence dans le vecteur indligneNA


enquete4_sans_NA<- enquete4[-indligneNA,] # Base de donnée sans les 69 valeurs manquantes
res.mca<- MCA(enquete4_sans_NA[,-c(1,2)],graph=FALSE,ncp = 36) # on effectue l'ACM

# valeur propres / variances
eig.val <- get_eigenvalue(res.mca) # 1/22= 0.04

# visualisation
fviz_screeplot (res.mca, addlabels = TRUE, ylim = c (0, 10)) # on garde entre 3 et 4 dimensions
fviz_eig(res.mca, choice = c("eigenvalue"), addlabels = TRUE)
# individus
ind <- get_mca_ind (res.mca)
fviz_mca_ind(res.mca,axes = c(1,2),
             label="none",
             repel = TRUE, 
             ggtheme = theme_minimal())
# Cos2 des individus
fviz_cos2(res.mca, choice = "ind", axes = 1:2, top = 30)
# Contribution des individus aux dimensions
fviz_contrib(res.mca, choice = "ind", axes = 1, top = 30)
fviz_contrib(res.mca, choice = "ind", axes = 2, top = 30)
fviz_contrib(res.mca, choice = "ind", axes = 3, top = 30)
fviz_contrib(res.mca, choice = "ind", axes = 1:2, top=30)
fviz_contrib(res.mca, choice = "ind", axes = 2:3, top = 30)
# graphques des variables
# Nous donne juste une catégorisation des variables
# Quelles questions sont plus influentes pour determiner la position de quelqun au long d'une dimension
# corrélation entre les dimensions et les questions : 
# Les questions plus centrées ne jouent aucun role dans les dimensions (ils peuvent repondre à n'importe quoi mais ça ne va pas changer)
fviz_mca_var (res.mca, choice = "mca.cor", axes = c(1,2), 
              repel = TRUE, 
              ggtheme = theme_minimal ())
fviz_mca_var(res.mca, col.var="contrib",axes = c(1,2),select.ind = list(contrib=10))

# Top 10 des variables actives avec le cos2 le plus elevé
fviz_mca_var (res.mca, select.var = list (cos2 = 10))

# Contribution totale aux dimensions 1 et 2
fviz_contrib(res.mca, choice = "var", axes =4, top = 30)

# Top 5 des categories de variables les plus contributives (meilleure visualisation)
# remarque : geom.ind or geom.var pour spécifier si on doit afficher les c(point, text or les deux)
fviz_mca_biplot (res.mca, select.ind = list (contrib = 402), axes = c(1,2), geom.ind = "point",
                 select.var = list (contrib = 20),
                 ggtheme = theme_minimal ())
# remarque : le rapprochement de deux individus dans le nouveau espace trouvé depdend cependant du cos^2 

fviz_ellipses(res.mca, c("Q3.4", "Q3.5","Q3.10"), select.ind = list(contrib=20),
              geom = "point")



############## ACM sans l'individu 224 : id : 11191749516 et sans données manquantes ##########################################

enquete4_sans_224<- enquete4_sans_NA[-185,]
res.mca.224<- MCA(enquete4_sans_224, graph = FALSE)

# visualisation
fviz_screeplot (res.mca.224, addlabels = TRUE, ylim = c (0, 10))
fviz_eig(res.mca.224, choice = c("eigenvalue"),addlabels = TRUE)
# individus
ind <- get_mca_ind (res.mca.224)
fviz_mca_ind(res.mca.224, col.ind = "cos2", label="none",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE, 
             ggtheme = theme_minimal())

# Variables
fviz_mca_var (res.mca.224, choice = "mca.cor", axes = c(1,2), 
              repel = TRUE, 
              ggtheme = theme_minimal ())
fviz_mca_var(res.mca.224,axes = c(1,4),col.var = "cos2")

fviz_mca_biplot (res.mca.224, select.ind = list (contrib = 402), axes = c(1,2),
                 select.var = list (contrib =30),
                 ggtheme = theme_minimal ())
###################### Données manquates traitement ############################################

install.packages("visdat")
library(visdat)
vis_miss(enquete3)
# autre façon 
install.packages("naniar")
library(naniar)
gg_miss_var(enquete4.new)



######### ACM sur la base de donnée 4 en laisant les données manquantes et en les imputants ##############
help(MCA)
tab = imputeMCA(enquete4, ncp=4)
res.mca.NA = MCA(tab$completeObs,graph=FALSE)

# valeur propres/ variances
eig.val <- get_eigenvalue(res.mca)

fviz_screeplot (res.mca.NA, addlabels = TRUE, ylim = c (0, 10))

# graphiques des individus 
fviz_mca_ind(res.mca.NA, col.ind = "cos2", label="none",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE, 
             ggtheme = theme_minimal())

# graphques des variables
var <- get_mca_var (res.mca.NA)
fviz_mca_var (res.mca.NA, choice = "mca.cor", axes = c(1,2),
              repel = TRUE, 
              ggtheme = theme_minimal ())

# Top 10 des variables actives avec le cos2 le plus elevé
fviz_mca_var (res.mca.NA, select.var = list (cos2 = 5))

# Contribution totale aux dimensions 1 et 2
fviz_contrib(res.mca.NA, choice = "var", axes = 2, top = 20)


# Top 5 des categories de variables les plus contributives (meilleure visualisation)
fviz_mca_biplot (res.mca.NA, select.ind = list ( contrib = 402 ),
                 select.var = list ( contrib = 30),
                 ggtheme = theme_minimal ())

################ ACM sans l'individu 224 avec données manquantes imputées #############################

help(imputeMCA)

enquete4_sans_224_NA<- enquete4[-223,]
is.na(enquete4_sans_224_NA)
tab.NA = imputeMCA(enquete4_sans_224_NA, ncp=4)
res.mca.NA.sans_224 = MCA(tab.NA$completeObs,graph=FALSE)

# visualisation
help(fviz_eig)
fviz_screeplot (res.mca.NA.sans_224, addlabels = TRUE, ylim = c (0, 10))

# individus
ind <- get_mca_ind (res.mca.NA.sans_224)
help(fviz_mca_ind)
fviz_mca_ind(res.mca.NA.sans_224, col.ind = "cos2", label= "none", axes = c(1,2),
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,
             ggtheme = theme_minimal())
help(MFA)
# graphques des variables
var <- get_mca_var (res.mca.224.NA)
help(fviz_mca_var)
fviz_mca_var (res.mca.NA.sans_224, choice = "mca.cor",
              repel = TRUE, 
              ggtheme = theme_minimal ())

# Contribution totale aux dimensions 1 et 2
fviz_contrib(res.mca.NA.sans_224, choice = "ind", axes = 3, top = 10)


# Top 5 des categories de variables les plus contributives (meilleure visualisation)
fviz_mca_biplot (res.mca.NA.sans_224, select.ind = list (contrib = 402), axes = c(2,3),
                 select.var = list (contrib = 20),
                 ggtheme = theme_minimal ())




############### Régroupement de certaines catégories mal représentées ############################

# questions : avis sur les pratiques collectifis
# nous allons avoir que trois catégories à savoir : Tout à fait d'accord, d'accord, pas d'accord
enquete3.new<- enquete3[,1:3]
enquete3.new$Q1.1.new <- fct_collapse(enquete3$Q1.1,
                                  "Pas d’accord" = c("Plutôt pas d’accord"),
                                  "D'accord" = c("Plutôt d’accord"))

enquete3.new$Q1.2.new <- fct_collapse(enquete3$Q1.2,
                                  "Pas d’accord" = c("Plutôt pas d’accord"),
                                  "D'accord" = c("Plutôt d’accord"))

# questions : avis sur les pratiques individuels

enquete3.new$Q2.1.new <- fct_collapse(enquete3$Q2.1,
                                  "Non" = c("Ne se prononce pas"))

enquete3.new$Q2.2.new <- fct_collapse(enquete3$Q2.2,
                                  "Non" = c("Ne se prononce pas"))

enquete3.new$Q2.3.new <- fct_collapse(enquete3$Q2.3,
                                  "Non" = c("Ne se prononce pas"))

enquete3.new$Q2.4.new <- fct_collapse(enquete3$Q2.4,
                                  "Pvicap" = c("Ne se prononce pas","Pas de voyages internationaux dans le cadre de mes activités professionnelles"),
                                  "La moitié environ (<20-60%)" = c("Pratiquement aucun en fait (<20%)","Quelques-uns (20-40%)","La moitié environ (40-60%)"))
enquete3.new$Q2.5<- enquete3[,10]
enquete3.new$Q2.6.new <- fct_collapse(enquete3$Q2.6,
                                  "Ne se prononce pas" = c("Je ne sais pas"))

enquete3.new$Q2.7.new <- fct_collapse(enquete3$Q2.7,
                                  "Oui, partiellement" = c("Oui, faiblement","Oui, moyennement"),
                                  "Non" = c("Ne se prononce pas"))

enquete3.new$Q2.8.new <- fct_collapse(enquete3$Q2.8,
                                  "Oui, partiellement" = c("Oui, faiblement","Oui, moyennement"),
                                  "Non" = c("Ne se prononce pas"))

enquete3.new$Q2.9.new <- fct_collapse(enquete3$Q2.9,
                                  "Oui, partiellement" = c("Oui, faiblement","Oui, moyennement"),
                                  "Non" = c("Ne se prononce pas"))

# Questions : Avis sur la position du cirad
enquete3.new$Q3.1.new <- fct_collapse(enquete3$Q3.1,
                                  "Pas d’accord" = c("Plutôt pas d’accord"),
                                  "d'accord" = c("Plutôt d’accord"))

enquete3.new$Q3.2.new <- fct_collapse(enquete3$Q3.2,
                                  "Pas d’accord" = c("Plutôt pas d’accord"),
                                  "d'accord" = c("Plutôt d’accord"))

enquete3.new$Q3.3.new <- fct_collapse(enquete3$Q3.3,
                                  "Pas d’accord" = c("Plutôt pas d’accord"),
                                  "d'accord" = c("Plutôt d’accord"))

enquete3.new$Q3.4.new <- fct_collapse(enquete3$Q3.4,
                                  "défavorable" = c("Tout à fait défavorable","Plutôt défavorable","Ne se prononce pas"))


enquete3.new$Q3.5.new <- fct_collapse(enquete3$Q3.5,
                                  "défavorable" = c("Tout à fait défavorable","Plutôt défavorable","Ne se prononce pas"))

enquete3.new$Q3.6.new <- fct_collapse(enquete3$Q3.6,
                                  "défavorable" = c("Tout à fait défavorable","Plutôt défavorable"))

enquete3.new$Q3.7.new <- fct_collapse(enquete3$Q3.7,
                                  "défavorable" = c("Tout à fait défavorable","Plutôt défavorable"))

enquete3.new$Q3.8.new <- fct_collapse(enquete3$Q3.8,
                                  "défavorable" = c("Tout à fait défavorable","Plutôt défavorable"))

enquete3.new$Q3.9.new <- fct_collapse(enquete3$Q3.9,
                                  "défavorable" = c("Tout à fait défavorable","Plutôt défavorable"))

enquete3.new$Q3.10.new <- fct_collapse(enquete3$Q3.10,
                                   "défavorable" = c("Tout à fait défavorable","Plutôt défavorable"))

enquete3.new$Q3.11.new <- fct_collapse(enquete3$Q3.11,
                                   "défavorable" = c("Tout à fait défavorable","Plutôt défavorable"))

enquete3.new$Q3.12.new <- fct_collapse(enquete3$Q3.12,
                                   "défavorable" = c("Tout à fait défavorable","Plutôt défavorable"))

enquete3.new$Q3.13.new <- fct_collapse(enquete3$Q3.13,
                                   "défavorable" = c("Tout à fait défavorable","Plutôt défavorable"))

# Questions : Avis socio-démographiques
enquete3.new$Q4.1.new<- enquete3[,28]
enquete3.new$Q4.2.new<- enquete3[,29]
enquete3.new$Q4.3.new <- fct_collapse(enquete3$Q4.3,
                                  "Chercheur" = c("cadre","collaborateur"), "Doctorant.e" = c("allocataire de recherche"),"Ingénieur,technicien-admin" = c("agent de maîtrise"))


enquete3.new$Q4.4.new<- fct_collapse(enquete3$Q4.4,
                                 "Chercheur" = c("Chercheuse","dans un autre institut","Maître de Conférences","chargé de recherche","Postdoc"),
                                 "Doctorant.e" = c("doctorant accueilli","doctorant","Doctorant","doctorante","Doctorante CIFRE","thesard","thésard"),
                                 "Stagiaire, VSC" = c("stagiaire","Stagiaire","VSC","Volontaire service civique"),
                                 "Ingénieur,technicien-admin" = c("Assistant ingénieur","CDD","ingénieur","ingenieur d'etude CDD","Ingénieur de recherche","ITA","technicien agricole et pépiniériste"))
enquete3.new$Q4.5.new<- enquete3[,32]
varName <- names(enquete3.new)[!names(enquete3.new)%in%c("Respondent.ID","Start.Date","End.Date")]
for (v in varName){
  enquete3.new[,v] <- as.character(enquete3.new[,v])
}


# Remplacer les information "Autre" dans Q4.3.new par celles du Q4.4.new
# et éliminer la colonne Q4.4
ind_ligne<- which(enquete3.new$Q4.4.new != "NA")
enquete3.new$Q4.3.new[ind_ligne]<- enquete3.new$Q4.4.new[ind_ligne]
enquete3.new<- enquete3.new[,-31]

# Traitement des réponses particulières en NA
enquete3.new$Q4.3.new[223]<- NA
enquete3.new$Q4.3.new[25]<- NA



#ggbivariate(data = enquete3.new, outcome = "Q4.1.new", explanatory = c("Q3.4.new","Q3.5.new","Q3.6.new","Q3.7.new","Q3.8.new","Q3.9.new","Q3.10.new"))


# nouvelle base de données avec les récatégorisation effectuées

enquete4.new<- enquete3.new[,-c(1,2,3)]




######### ACM sur la base de données enquete 4 en elevant les données manquantes #####################


# traitement de données manquantes
is.na(enquete4.new) # nombre de données manquantes
# on enlève de la base les individus ayant au moins un NA
indligneNA.new <- which(is.na(enquete4.new),arr.ind = TRUE)[,1] # indice des lignes dans la colonne 1
length(unique(indligneNA.new)) # 69 individus ont au moins une donnée manquante

table(indligneNA.new) # nombre d'occurence dans le vecteur indligneNA.new

enquete4_sans_NA.new<- enquete4.new[-indligneNA.new,]
res.mca.new<- MCA(enquete4_sans_NA.new, graph=FALSE,
                  quali.sup = 25:28)

# valeur propres / variances
eig.val.new <- get_eigenvalue(res.mca.new) # 1/22= 0.04

# visualisation
help(fviz_eig)
fviz_screeplot (res.mca.new, addlabels = TRUE, ylim = c (0, 10))
fviz_eig(res.mca.new, choice = c("eigenvalue"), addlabels = TRUE) # on garde deux axes principaux
# individus
ind <- get_mca_ind (res.mca.new)

fviz_mca_ind(res.mca.new, col.ind = "cos2",

             label = "none", axes = c(1,2),
             repel = TRUE, 
             ggtheme = theme_minimal())

fviz_mca_ind(res.mca.new,
             habillage = "Q4.1.new", label = "none",
             repel = TRUE, ellipse.type = "confidence", addEllipses = TRUE,
             ggtheme = theme_minimal())
# Cos2 des individus
fviz_cos2(res.mca.new, choice = "ind", axes = 1:2, top = 20)
# Contribution des individus aux dimensions
fviz_contrib(res.mca.new, choice = "ind", axes = 1:2, top = 20)

# graphques des variables
fviz_mca_var (res.mca.new, choice = "mca.cor", axes = c(2,3),
              repel = TRUE, 
              ggtheme = theme_minimal ())


# Top 10 des variables actives avec le cos2 le plus elevé
fviz_mca_var (res.mca.new, select.var = list (cos2 = 20),axes = c(3,4))

# Contribution totale aux dimensions 1 et 2
fviz_contrib(res.mca.new, choice = "var", axes =3, top = 25)

# Top 5 des categories de variables les plus contributives (meilleure visualisation)
fviz_mca_biplot (res.mca.new, select.ind = list (contrib = 402), axes = c(2,3), geom.ind = "point",
                 select.var = list (contrib = 20),
                 ggtheme = theme_minimal ())
# remarque : le rapprochement de deux individus dans le nouveau espace trouvé depdend cependant du cos^2 



############## ACM sans l'individu 224 : id : 11191749516 et sans données manquantes news ##########################################

enquete4_sans_224.new<- enquete4_sans_NA.new[-185,]
res.mca.224.new<- MCA(enquete4_sans_224.new, graph = FALSE)

# visualisation
help(fviz_eig)
fviz_screeplot (res.mca.224.new, addlabels = TRUE, ylim = c (0, 10))
# on garde aussi deux axes 
# individus
ind.new <- get_mca_ind (res.mca.224.new)
fviz_mca_ind(res.mca.224.new, col.ind = "cos2", label = "none",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE, 
             ggtheme = theme_minimal())

# Top 5 des categories de variables les plus contributives (meilleure visualisation)
fviz_mca_biplot (res.mca.224.new, select.ind = list (contrib = 402), axes = c(2,3),
                 select.var = list (contrib = 20),
                 ggtheme = theme_minimal ())

fviz_contrib(res.mca.224.new,choice = "ind",axes = 3,top = 30)

######### ACM sur la base de donnée 4 new en laisant les données manquantes et en les imputants ##############


help(imputeMCA)
tab.comp.new = imputeMCA(enquete4.new, ncp=4)
res.mca.NA.new = MCA(tab.comp.new$completeObs,graph=FALSE, quali.sup = 25:28)

# valeur propres/ variances
#eig.val <- get_eigenvalue(res.mca)

fviz_screeplot (res.mca.NA.new, addlabels = TRUE, ylim = c (0, 10))
# 2 aussi axes à garder 

# graphiques des individus 
fviz_mca_ind(res.mca.NA.new, col.ind = "cos2", label="none",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE, 
             ggtheme = theme_minimal())

# graphques des variables
var <- get_mca_var (res.mca.NA)
fviz_mca_var (res.mca.NA.new, choice = "mca.cor",
              repel = TRUE, 
              ggtheme = theme_minimal ())


# Top 5 des categories de variables les plus contributives (meilleure visualisation)
fviz_mca_biplot (res.mca.NA.new, select.ind = list (contrib = 402), axes = c(2,3),
                 select.var = list (contrib = 20),
                 ggtheme = theme_minimal ())


################ ACM sans l'individu 224 avec données manquantes imputées #############################

# On enlève surtout l'individu extrême mais aussi tous les individus qui ont au moins
# 19 valeurs manquantes sur 28 variables avnt d'effectuer l'imputation

enquete4_sans_224_NA.new<- enquete4.new[-c(223,17,26,31,33,37,38,47,54,
                                           58,97,108,118,125,126,133,157,160,193
                                           ,201,203,204,209,228,234,244,245,265,
                                           267,290,330,333,337,352,354,378,379,381),]
tab.comp.NA.new = imputeMCA(enquete4_sans_224_NA.new)
res.mca.224.NA.new<- MCA(tab.comp.NA.new$completeObs, graph = FALSE, quali.sup = 25:28)

# Statistique descriptives


ggbivariate(data = tab.comp.NA.new$completeObs, outcome = "Q4.1.new", explanatory = c("Q3.4.new","Q3.5.new","Q3.6.new","Q3.7.new","Q3.8.new","Q3.9.new","Q3.10.new"))
ggtable(data = tab.comp.NA.new$completeObs,  columnsX =   "Q4.1.new",columnsY = c("Q1.1.new","Q1.2.new"))

# visualisation

fviz_screeplot (res.mca.224.NA.new, addlabels = TRUE, ylim = c (0, 10))

# individus

fviz_mca_ind(res.mca.224.NA.new, axes = c(2,3),
             habillage = "Q2.5.new", label="none",
             repel = TRUE, addEllipses = TRUE, ellipse.type="confidence",
             ggtheme = theme_minimal())

# graphques des variables
var <- get_mca_var (res.mca.224.NA)
fviz_mca_var (res.mca.224.NA.new, choice = "mca.cor",
              repel = TRUE, axes = c(2,3),
              ggtheme = theme_minimal ())

# Contribution totale aux dimensions - variables
fviz_contrib(res.mca.224.NA.new, choice = "var", axes = 1, top = 50)
fviz_contrib(res.mca.224.NA.new, choice = "var", axes = 2, top = 50)
fviz_contrib(res.mca.224.NA.new, choice = "var", axes = 3, top = 50)

fviz_contrib(res.mca.224.NA.new, choice = "var", axes = 4, top = 50)

# contribution totale aux dimensions - individus
fviz_cos2(res.mca.224.NA.new, choice = "ind", axes = 1, top = 50)
fviz_cos2(res.mca.224.NA.new, choice = "ind", axes = 2, top = 50)
fviz_cos2(res.mca.224.NA.new, choice = "ind", axes = 3, top = 50)

fviz_cos2(res.mca.224.NA.new, choice = "ind", axes = 4, top = 50)

# Biplot individus-variables
fviz_mca_biplot (res.mca.224.NA.new, select.ind = list (contrib = 402), axes = c(2,3),
                 select.var = list (contrib = 100),
                 ggtheme = theme_minimal ())




############ ACP ####################################
enquete3.mixte<- enquete3.new[,1:3]

# Avis sur les pratiques colletif
enquete3.mixte$Q1.1.quant<- fct_recode(enquete3.new$Q1.1.new,
                                       "-1" = "Ne se prononce pas",
                                       "0" = "Pas d’accord",
                                       "1" = "D'accord",
                                       "2" = "Tout à fait d’accord")

enquete3.mixte$Q1.2.quant<- fct_recode(enquete3.new$Q1.2.new,
                                       "-1" = "Ne se prononce pas",
                                       "0" = "Pas d’accord",
                                       "1" = "D'accord",
                                       "2" = "Tout à fait d’accord")
# Avis sur les pratiques individiduels

enquete3.mixte$Q2.1.quant<- fct_recode(enquete3.new$Q2.1.new,
                                       "-1" = "Pas concerné.e",
                                       "0" = "Non",
                                       "1" = "Oui, légèrement",
                                       "2" = "Oui, significativement")

enquete3.mixte$Q2.2.quant<- fct_recode(enquete3.new$Q2.2.new,
                                       "0" = "Non",
                                       "1" = "Oui, légèrement",
                                       "2" = "Oui, significativement")

enquete3.mixte$Q2.3.quant<- fct_recode(enquete3.new$Q2.3.new,
                                       "0" = "Non",
                                       "1" = "Oui")

enquete3.mixte$Q2.4.quant<- fct_recode(enquete3.new$Q2.4.new,
                                       "-1" = "Pvicap",
                                       "0" = "La moitié environ (<20-60%)",
                                       "1" = "Une grande partie (60-80%)",
                                       "2" = "Tous ou presque (>80%)")

enquete3.mixte$Q2.5.quant<- fct_recode(enquete3.new$Q2.5,
                                       "0" = "Ne se prononce pas",
                                       "1" = "Plutôt non",
                                       "2" = "Plutôt oui")

enquete3.mixte$Q2.6.quant<- fct_recode(enquete3.new$Q2.6.new,
                                       
                                       "-1" = "Ne se prononce pas",
                                       "0" = "Non",
                                       "1" = "Peut-être/Ça dépend",
                                       "2" = "Oui")

enquete3.mixte$Q2.7.quant<- fct_recode(enquete3.new$Q2.7.new,
                               
                                       "0" = "Non",
                                       "1" = "Oui, partiellement",
                                       "2" = "Oui, fortement")


enquete3.mixte$Q2.8.quant<- fct_recode(enquete3.new$Q2.8.new,
                                       
                                       "0" = "Non",
                                       "1" = "Oui, partiellement",
                                       "2" = "Oui, fortement")


enquete3.mixte$Q2.9.quant<- fct_recode(enquete3.new$Q2.9.new,
                                       
                                       "0" = "Non",
                                       "1" = "Oui, partiellement",
                                       "2" = "Oui, fortement")

# Avis sur la position du CIRAD


enquete3.mixte$Q3.1.quant<- fct_recode(enquete3.new$Q3.1.new,
                                       "-1" = "Ne se prononce pas",
                                       "0" = "Pas d’accord",
                                       "1" = "d'accord",
                                       "2" = "Tout à fait d’accord")

enquete3.mixte$Q3.2.quant<- fct_recode(enquete3.new$Q3.2.new,
                                       "-1" = "Ne se prononce pas",
                                       "0" = "Pas d’accord",
                                       "1" = "d'accord",
                                       "2" = "Tout à fait d’accord")

enquete3.mixte$Q3.3.quant<- fct_recode(enquete3.new$Q3.3.new,
                                       "-1" = "Ne se prononce pas",
                                       "0" = "Pas d’accord",
                                       "1" = "d'accord",
                                       "2" = "Tout à fait d’accord")

enquete3.mixte$Q3.4.quant<- fct_recode(enquete3.new$Q3.4.new,
                                       
                                       "0" = "défavorable",
                                       "1" = "Plutôt favorable",
                                       "2" = "Tout à fait favorable")

enquete3.mixte$Q3.5.quant<- fct_recode(enquete3.new$Q3.5.new,
                                      
                                       "0" = "défavorable",
                                       "1" = "Plutôt favorable",
                                       "2" = "Tout à fait favorable")

enquete3.mixte$Q3.6.quant<- fct_recode(enquete3.new$Q3.6.new,
                                       "-1" = "Ne se prononce pas",
                                       "0" = "défavorable",
                                       "1" = "Plutôt favorable",
                                       "2" = "Tout à fait favorable")

enquete3.mixte$Q3.7.quant<- fct_recode(enquete3.new$Q3.7.new,
                                       "-1" = "Ne se prononce pas",
                                       "0" = "défavorable",
                                       "1" = "Plutôt favorable",
                                       "2" = "Tout à fait favorable")

enquete3.mixte$Q3.8.quant<- fct_recode(enquete3.new$Q3.8.new,
                                       "-1" = "Ne se prononce pas",
                                       "0" = "défavorable",
                                       "1" = "Plutôt favorable",
                                       "2" = "Tout à fait favorable")

enquete3.mixte$Q3.9.quant<- fct_recode(enquete3.new$Q3.9.new,
                                       "-1" = "Ne se prononce pas",
                                       "0" = "défavorable",
                                       "1" = "Plutôt favorable",
                                       "2" = "Tout à fait favorable")

enquete3.mixte$Q3.10.quant<- fct_recode(enquete3.new$Q3.10.new,
                                       "-1" = "Ne se prononce pas",
                                       "0" = "défavorable",
                                       "1" = "Plutôt favorable",
                                       "2" = "Tout à fait favorable")

enquete3.mixte$Q3.11.quant<- fct_recode(enquete3.new$Q3.11.new,
                                       "-1" = "Ne se prononce pas",
                                       "0" = "défavorable",
                                       "1" = "Plutôt favorable",
                                       "2" = "Tout à fait favorable")

enquete3.mixte$Q3.12.quant<- fct_recode(enquete3.new$Q3.12.new,
                                       "-1" = "Ne se prononce pas",
                                       "0" = "défavorable",
                                       "1" = "Plutôt favorable",
                                       "2" = "Tout à fait favorable")

enquete3.mixte$Q3.13.quant<- fct_recode(enquete3.new$Q3.13.new,
                                       "-1" = "Ne se prononce pas",
                                       "0" = "défavorable",
                                       "1" = "Plutôt favorable",
                                       "2" = "Tout à fait favorable")

enquete3.mixte$Q4.1.quant<- enquete3.new[,28]
enquete3.mixte$Q4.2.quant<- enquete3.new[,29]
enquete3.mixte$Q4.3.quant<- enquete3.new[,30]
enquete3.mixte$Q4.5.quant<- enquete3.new[,31]

# conversion des variables en numérique
varName.mixte <- names(enquete3.mixte)[!names(enquete3.mixte)%in%c("Respondent.ID","Start.Date","End.Date","Q4.1.quant","Q4.2.quant","Q4.3.quant","Q4.5.quant")]
for (v in varName.mixte){
  enquete3.mixte[,v] <- as.numeric(as.character(enquete3.mixte[,v]))
}

enquete4.quant<- enquete3.mixte[-223,-c(1,2,3)]

#### ACP sur la base de donnée quantitative sans données manquantes ##################
# on enlève de la base les individus ayant au moins un NA
indligneNA.mixte <- which(is.na(enquete4.quant),arr.ind = TRUE)[,1] # indice des lignes dans la colonne 1
length(unique(indligneNA.mixte)) # 69 individus ont au moins une donnée manquante
enquete4.quant_sans_NA<- enquete4.quant[-indligneNA.mixte,]

res.pca.sans_NA<- PCA(enquete4.quant_sans_NA,graph=FALSE, quali.sup = 25:28)

# Valeurs propres

eig.val <- get_eigenvalue(res.pca.sans_NA)
eig.val

fviz_eig(res.pca.sans_NA, addlabels = TRUE, ylim = c(0, 50))

# Variables

# Colorer en fonction du cos2: qualité de représentation
fviz_pca_var(res.pca.sans_NA, col.var = "cos2", axes = c(2,3),
             repel = TRUE # Évite le chevauchement de texte
)


# Contributions des variables à PC1
fviz_contrib(res.pca.sans_NA, choice = "var", axes = 2, top = 20)
# Contributions des variables à PC2
fviz_contrib(res.pca.sans_NA, choice = "var", axes = c(1,2), top = 20)


# Individus
fviz_contrib(res.pca.sans_NA, choice = "ind", axes = 1, top=50)

fviz_pca_ind (res.pca.sans_NA, label="none",
              repel = TRUE # Évite le chevauchement de texte
)
# biplot
fviz_pca_biplot(res.pca.sans_NA, repel = TRUE, geom.ind = "point", axes = c(1,2),
                col.var = "#2E9FDF", # Couleur des variables
                col.ind = "#696969"  # Couleur des individues
)

fviz_pca_biplot (res.pca.sans_NA, select.ind = list (contrib = 402), axes = c(1,2),
                 select.var = list (contrib = 30), geom.ind = "point",
                 ggtheme = theme_minimal())



################## MFA ##################################

enquete4.MFA<- enquete3.new[,-c(1,2,3)]

######### avec imputation des données manquantes ##########################
enquete4_sans_224_NA.MFA<- enquete4.MFA[-c(223,17,26,31,33,37,38,47,54,
                                           58,97,108,118,125,126,133,157,160,193
                                           ,201,203,204,209,228,234,244,245,265,
                                           267,290,330,333,337,352,354,378,379,381),]



tab.comp.NA.MFA = imputeMCA(enquete4_sans_224_NA.MFA)

res.NA.MFA<- MFA(tab.comp.NA.MFA$completeObs, group = c(2,9,13,4),
                 name.group = c("Pratiques collectives","Pratiques individuelles",
                                "Position CIRAD","socio-démographique"),
                 type=c("n","n","n","n"),num.group.sup = 4,graph=FALSE)


# Valeur propres

fviz_screeplot(res.NA.MFA,addlabels=TRUE)
fviz_eig(res.NA.MFA, choice = c("eigenvalue"), addlabels = TRUE)

# Graphiques des variables en fonction des groupes
group<- get_mfa_var(res.NA.MFA, "group");group
head(group$coord)
head(group$contrib)

# corrélation entre les groupes et le plan (1,2)
fviz_mfa_var(res.NA.MFA, "group",axes = c(1,2))

# corrélation entre les groupes et le plan (2,3)
fviz_mfa_var(res.NA.MFA, "group",axes = c(2,3))

# corrélation entre les groupes et le plan (3,4)
fviz_mfa_var(res.NA.MFA, "group",axes = c(3,4), invisible="quali.var")

# corrélation entre les groupes et le plan (4,5)
fviz_mfa_var(res.NA.MFA, "group",axes = c(4,5))

# Individus


fviz_contrib(res.NA.MFA, "ind", axes = 2, top = 20)


fviz_cos2(res.NA.MFA, "ind", axes = 1, top = 60)

# nuage de point : Individus du point de vu de l'ensemble des groupes

fviz_mfa_ind(res.NA.MFA, col.ind =  "cos2", axes = c(3,4), invisible = "quali.var",
             geom.ind="point")
# en fonction des variables à l'étude

fviz_mfa_ind(res.NA.MFA, 
             habillage = "Q4.1.new", # color by groups 
             axes = c(2,3), 
             addEllipses = TRUE, ellipse.type = "confidence", invisible = "quali.var"
             
) 

# Individus partiels

fviz_mfa_ind (res.NA.MFA,repel = TRUE, label="none", partial= c(1,2,3,4,5))
fviz_mfa_quali_biplot(res.NA.MFA, geom.ind="point",
                      select.ind = list (contrib = 402), axes = c(1,2),
                      select.var = list (contrib = 70))


fviz_mfa_ind (res.NA.MFA, repel = TRUE, partial = c ("171","203","369", "62","388"))




#### on enlève de la base les individus ayant au moins un NA ######################
indligneNA.mfa <- which(is.na(enquete4.MFA),arr.ind = TRUE)[,1] # indice des lignes dans la colonne 1
length(unique(indligneNA.mfa)) # 74 individus ont au moins une donnée manquante
enquete4.mfa_sans_NA<- enquete4.MFA[-indligneNA.mfa,]


res.MFA<- MFA(enquete4.mfa_sans_NA, group = c(2,9,13,4),
              name.group = c("Pratiques collectifs","Pratiques individuels",
                             "Position CIRAD","socio-démographique"),
              type=c("n","n","n","n"),num.group.sup = 4,graph=FALSE)

# Valeur propres

fviz_screeplot(res.NA.MFA,addlabels=TRUE)
fviz_eig(res.MFA, choice = c("eigenvalue"), addlabels = TRUE)

# Graphiques des variables 
group<- get_mfa_var(res.MFA, "group");group
head(group$coord)

fviz_mfa_var(res.MFA, "group",axes = c(3,4))


# Individus


fviz_contrib(res.MFA, "ind", axes = 2, top = 20)


fviz_cos2(res.MFA, "ind", axes = 1, top = 60)


# en fonction des variables à l'étude
help("MFA")
fviz_mfa_ind(res.MFA, 
             habillage = "Q2.3.new", # color by groups 
              axes = c(4,5), 
             addEllipses = TRUE, ellipse.type = "confidence", invisible = "quali.var"
             
) 
#biplot
fviz_mfa_quali_biplot(res.MFA, repel = TRUE,   ggtheme = theme_minimal (), geom.ind="point",
                      select.ind = list (contrib = 402), axes = c(4,5),
                      select.var = list (contrib = 70))



# Individus partiels

fviz_mfa_ind (res.NA.MFA,repel = TRUE, label="none", partial= c(1,2,3,4,5))
fviz_mfa_quali_biplot(res.NA.MFA, repel = TRUE,   ggtheme = theme_minimal (), geom.ind="point",
                      select.ind = list (contrib = 402), axes = c(4,5),
                      select.var = list (contrib = 70))


fviz_mfa_ind (res.NA.MFA, repel = TRUE, partial = c ("171","203","369", "62","388"))






################### Méthodes de classifications non supervisées ##################################

library(cluster)
ppp.diss <- daisy(tab.comp.NA.new$completeObs, metric = "gower",type = list(ordratio = 1:24))
#ppp.diss<- as.matrix(ppp.diss)

# CAH à tester 

# en utilisant la distance de ward
enquete_cah<- hclust(ppp.diss, method = "ward.D2")
plot(enquete_cah)
clusters <- rect.hclust(enquete_cah, k=5)

print(sort(cutree(enquete_cah,k=5)))


### CAH en utilsant les résultats de MCA et MFA  ACP ##############

clu_acm<- HCPC(res.mca.224.NA.new, graph= FALSE) # ACM

clu_acp<- HCPC(res.pca.sans_NA, graph=FALSE) # ACP

clu_mfa_imputation<- HCPC(res.NA.MFA, graph=FALSE) # MFA avec imputation

clu_mfa<- HCPC(res.MFA, graph=FALSE) # MFA en enlèvant les 74 valeurs manquantes

# nombre d'individus dans chaque cluster
sum(clu_mfa$data.clust$clust==1)
sum(clu_mfa_imputation$data.clust$clust==1)

sum(clu_mfa_imputation$data.clust$clust==2)
sum(clu_mfa$data.clust$clust==2)

sum(clu_mfa_imputation$data.clust$clust==3)
sum(clu_mfa$data.clust$clust==3)

sum(clu_mfa_imputation$data.clust$clust==4)
sum(clu_mfa$data.clust$clust==4)

sum(clu_mfa_imputation$data.clust$clust==5)



# Dendrogramme
# MFA avec 364 observations
fviz_dend(clu_mfa_imputation, show_labels = FALSE) # 5 classes mais on peut chosir 4 classes

#MFA avec 328 observations
fviz_dend(clu_mfa, show_labels = FALSE) # 4 classes 

# ACM avec 364 observations
fviz_dend(clu_acm, show_labels = FALSE) # 3 classes avec l'ACM

# ACP avec 328 observations
fviz_dend(clu_acp, show_labels = FALSE) # 3 classes avec l'ACP

# Individus
options(ggrepel.max.overlaps = Inf)
fviz_cluster(clu_mfa_imputation, geom = "text", main = "Factor map", repel=TRUE)

fviz_cluster(clu_mfa, geom = "text", main = "Factor map", repel=FALSE) + scale_color_manual(values=c("green", "violet", "blue","red")) 

fviz_cluster(clu_acm, geom = "text", main = "Factor map", repel=TRUE)

fviz_cluster(clu_acp, geom="text", main = "Factor map", repel=TRUE)


# graphe nuage de points en fonction des clusters ACP

donne<- clu_acp$data.clust
res0<- PCA(donne,graph=FALSE, quali.sup = 25:29)

fviz_pca_ind (res0, habillage="clust",axes=c(1,2),addEllipses=TRUE,
              ellipse.type="confidence", invisible="quali.var"
)

# graphe nuage de points en fonction des clusters obtenus par ACM
donne1<- clu_acm$data.clust
res1<- MCA(donne1, graph = FALSE, quali.sup = 25:29)

fviz_mca_ind(res1, 
             habillage = "clust", # color by groups 
             axes = c(2,3), 
             addEllipses = TRUE, ellipse.type = "confidence", invisible = "quali.var"
             
) 
# graphe nuage de points en fontion des clusters obtenus par MCA

# clustering sur les questions du groupe position du CIRAD

clus_Cirad<- HCPC(res.NA.MFA$separate.analyses$`Pratiques collectives`, graph = FALSE)
fviz_dend(clus_Cirad, show_labels = FALSE)
donne4<- clus_Cirad$data.clust
res4<- MCA(donne4, graph = FALSE)
fviz_mca_ind(res4, 
             habillage = "clust", # color by groups 
             axes = c(1,2), 
             addEllipses = TRUE, ellipse.type = "confidence", invisible = "quali.var"
) 

donne2<- clu_mfa$data.clust

write.table(donne2[,29], 
            "/Users/abdoulayebaradji/Desktop/par_question/clus.txt", sep= " ", eol="\n")


res2<- MFA(donne2, group = c(2,9,13,5),
                name.group = c("Pratiques collectifs","Pratiques individuels",
                               "Position CIRAD","socio-démographique"),
                type=c("n","n","n","n"),num.group.sup = 4,graph=FALSE)
fviz_mfa_ind(res2, 
             habillage = "clust", # color by groups 
             axes = c(1,2), 
             addEllipses = TRUE, ellipse.type = "confidence", invisible = "quali.var"
)





########## Analyse des questions ouvertes ##################################

# On exporte pour chaque colonne, les réponses aux questions ouvertes

enquete_ouverte<- enquete1[,c(1,5,7,9,11,13,15,17,19,21,23,25,27,29,31,42)]
names(enquete_ouverte) <- c("Respondent ID",                                                                                                                        
 "c_Q1.1"              ,                                                                                                                   
 "c_Q1.2"            ,                                                                                                                    
 "c_Q2.1"   ,                                                                                                                     
 "c_Q2.2"  ,                                                                                                                      
 "c_Q2.3"        ,                                                                                                                       
 "c_Q2.4"    ,                                                                                                                        
 "c_Q2.5"      ,                                                                                                                         
 "c_Q2.6" ,                                                                                                                          
 "c_Q2.7",
 "c_Q2.8"   ,                                                                                                                            
 "c_Q2.9"  ,                                                                                                                             
"c_Q3.1",
"c_Q3.2",
"c_Q3.3",
 "Avez-vous des commentaires ou des suggestions en lien avec les actions possibles pour réduire l’empreinte environnementale du CIRAD ?"
)
l<- c(which(enquete_ouverte[2]!="NA"), which(enquete_ouverte[3]!="NA"))
l<- unique(l)

write.table(enquete_ouverte[,1], 
            "/Users/abdoulayebaradji/Desktop/par_question/id.txt", sep= " ", eol="\n")
write.table(enquete_ouverte[which(enquete_ouverte[2]!="NA"),1], 
            "/Users/abdoulayebaradji/Desktop/par_question/id_Q1.1.txt", sep= " ", eol="\n")

write.table(enquete_ouverte[l,2:3], "/Users/abdoulayebaradji/Desktop/par_question/Q1.union.txt", sep=" ", eol="\n")

write.table(enquete_ouverte[which(enquete_ouverte[2] !="NA"),2]
            ,"/Users/abdoulayebaradji/Desktop/par_question/Q1.1.txt", sep = " ", eol = "\n"
)

write.table(enquete_ouverte[which(enquete_ouverte[3]!="NA"),1], 
            "/Users/abdoulayebaradji/Desktop/par_question/id_Q1.2.txt", sep= " ", eol="\n")

write.table(enquete_ouverte[-1,3][which(enquete_ouverte[3] !="NA")]
            ,"/Users/abdoulayebaradji/Desktop/par_question/Q1.2.txt", sep = " ", eol = "\n"
)

write.table(enquete_ouverte[which(enquete_ouverte[4]!="NA"),1], 
            "/Users/abdoulayebaradji/Desktop/par_question/id_Q2.1.txt", sep= " ", eol="\n")

write.table(enquete_ouverte[-1,4][which(enquete_ouverte[4] !="NA")]
            ,"/Users/abdoulayebaradji/Desktop/par_question/Q2.1.txt", sep = " ", eol = "\n"
)

write.table(enquete_ouverte[which(enquete_ouverte[5]!="NA"),1], 
            "/Users/abdoulayebaradji/Desktop/par_question/id_Q2.2.txt", sep= " ", eol="\n")


write.table(enquete_ouverte[-1,5][which(enquete_ouverte[5] !="NA")]
            ,"/Users/abdoulayebaradji/Desktop/par_question/Q2.2.txt", sep = " ", eol = "\n"
)

write.table(enquete_ouverte[which(enquete_ouverte[6]!="NA"),1], 
            "/Users/abdoulayebaradji/Desktop/par_question/id_Q2.3..txt", sep= " ", eol="\n")


write.table(enquete_ouverte[-1,6][which(enquete_ouverte[6] !="NA")]
            ,"/Users/abdoulayebaradji/Desktop/par_question/Q2.3.txt", sep = " ", eol = "\n"
)

write.table(enquete_ouverte[which(enquete_ouverte[7]!="NA"),1], 
            "/Users/abdoulayebaradji/Desktop/par_question/id_Q2.4.txt", sep= " ", eol="\n")

write.table(enquete_ouverte[-1,7][which(enquete_ouverte[7] !="NA")]
            ,"/Users/abdoulayebaradji/Desktop/par_question/Q2.4.txt", sep = " ", eol = "\n"
)

write.table(enquete_ouverte[which(enquete_ouverte[8]!="NA"),1], 
            "/Users/abdoulayebaradji/Desktop/par_question/id_Q2.5.txt", sep= " ", eol="\n")

write.table(enquete_ouverte[-1,8][which(enquete_ouverte[8] !="NA")]
            ,"/Users/abdoulayebaradji/Desktop/par_question/Q2.5.txt", sep = " ", eol = "\n"
)

write.table(enquete_ouverte[which(enquete_ouverte[9]!="NA"),1], 
            "/Users/abdoulayebaradji/Desktop/par_question/id_Q2.6.txt", sep= " ", eol="\n")

write.table(enquete_ouverte[-1,9][which(enquete_ouverte[9] !="NA")]
            ,"/Users/abdoulayebaradji/Desktop/par_question/Q2.6.txt", sep = " ", eol = "\n"
)

write.table(enquete_ouverte[which(enquete_ouverte[10]!="NA"),1], 
            "/Users/abdoulayebaradji/Desktop/par_question/id_Q2.7.txt", sep= " ", eol="\n")

write.table(enquete_ouverte[-1,10][which(enquete_ouverte[10] !="NA")]
            ,"/Users/abdoulayebaradji/Desktop/par_question/Q2.7.txt", sep = " ", eol = "\n"
)

write.table(enquete_ouverte[which(enquete_ouverte[11]!="NA"),1], 
            "/Users/abdoulayebaradji/Desktop/par_question/id_Q2.8.txt", sep= " ", eol="\n")

write.table(enquete_ouverte[-1,11][which(enquete_ouverte[11] !="NA")]
            ,"/Users/abdoulayebaradji/Desktop/par_question/Q2.8.txt", sep = " ", eol = "\n"
)

write.table(enquete_ouverte[which(enquete_ouverte[12]!="NA"),1], 
            "/Users/abdoulayebaradji/Desktop/par_question/id_Q2.9.txt", sep= " ", eol="\n")

write.table(enquete_ouverte[-1,12][which(enquete_ouverte[12] !="NA")]
            ,"/Users/abdoulayebaradji/Desktop/par_question/Q2.9.txt", sep = " ", eol = "\n"
)

write.table(enquete_ouverte[which(enquete_ouverte[13]!="NA"),1], 
            "/Users/abdoulayebaradji/Desktop/par_question/id_Q3.1.txt", sep= " ", eol="\n")


write.table(enquete_ouverte[-1,13][which(enquete_ouverte[13] !="NA")]
            ,"/Users/abdoulayebaradji/Desktop/par_question/Q3.1.txt", sep = " ", eol = "\n"
)

write.table(enquete_ouverte[which(enquete_ouverte[14]!="NA"),1], 
            "/Users/abdoulayebaradji/Desktop/par_question/id_Q3.2.txt", sep= " ", eol="\n")

write.table(enquete_ouverte[which(enquete_ouverte[14] !="NA"),14]
            ,"/Users/abdoulayebaradji/Desktop/par_question/Q3.2.txt", sep = " ", eol = "\n"
)

write.table(enquete_ouverte[which(enquete_ouverte[15]!="NA"),1], 
            "/Users/abdoulayebaradji/Desktop/par_question/id_Q3.3.txt", sep= " ", eol="\n")

write.table(enquete_ouverte[-1,15][which(enquete_ouverte[15] !="NA"),]
            ,"/Users/abdoulayebaradji/Desktop/par_question/Q3.3.txt", sep = " ", eol = "\n"
)

write.table(enquete_ouverte[-1,16][which(enquete_ouverte[16] !="NA")]
            ,"/Users/abdoulayebaradji/Desktop/par_question/Q.txt", sep = " ", eol = "\n"
)

### par individus 

# cluster 4
write.table(enquete_ouverte[c(79,108,213,172,112),], na="NA", sep = "\n"
            ,"/Users/abdoulayebaradji/Desktop/par_individu/cluster4/clu.txt"
)

# cluster 3
write.table(enquete_ouverte[c(306,343,99,178,138,182,128,289,200,276,125,46,29,
                      188,197,295,325,242,313,299,231,184,150,223,312,176),], na="NA"
            ,"/Users/abdoulayebaradji/Desktop/par_individu/cluster3/clu_3.txt", sep = "\n"
)

# cluster 2


# cluster 1
write.table(enquete_ouverte[18,], na=" "
            ,"/Users/abdoulayebaradji/Desktop/par_individu/ind_18.txt", sep = " "
)




########### text mining sur les commentaires du manifeste ###################
install.packages("tm")
library(tm)
install.packages("readtext")
library(readtext)
library("SnowballC")
install.packages("wordcloud")
library("wordcloud")
library("RColorBrewer")

# importation donne
text<- readLines("data/commentaires.txt")

doc_ids <- c(1:82)

df <- data.frame(doc_id = doc_ids, text = text, stringsAsFactors = FALSE)
df_corpus <- Corpus(DataframeSource(df))

# Convertir le texte en minuscule
df_corpus <- tm_map(df_corpus, content_transformer(tolower))
# Supprimer les nombres
df_corpus <- tm_map(df_corpus, removeNumbers)
# Supprimer les ponctuations
df_corpus <- tm_map(df_corpus, removePunctuation)

# Supprimer les espaces vides supplémentaires
df_corpus <- tm_map(df_corpus, stripWhitespace)

# Supprimer les mots vides français
df_corpus <- tm_map(df_corpus, removeWords, stopwords("french"))
# Supprimer votre propre liste de mots non désirés
#docs <- tm_map(docs, removeWords, c("blabla1", "blabla2")) 

# Text stemming
#df_corpus <- tm_map(df_corpus, stemDocument)

dtm <- TermDocumentMatrix(df_corpus)
m <- as.matrix(dtm)
m<- t(m)
v <- sort(rowSums(m),decreasing=TRUE)
d <- data.frame(word = names(v),freq=v)
head(d, 10)

# nuage de mot

set.seed(1234)
wordcloud(words = d$word, freq = d$freq, min.freq = 5,
          max.words=200, random.order=FALSE, rot.per=0.35, 
          colors=brewer.pal(8, "Dark2"))
# Mots les plus fréquents
findFreqTerms(dtm, lowfreq = 20)

######## text mining en fonction des commentaires ou suggestion pour réduire l'empreinte environnementale du Cirad #################

install.packages("quanteda.textmodels")
library(quanteda.textmodels)
#library(quanteda_options(language = "fr"))
library(textclean)
library(stringr)

# importation donne
text_Q<- readLines("data/Q.txt")

doc_ids_Q <- c(1:34)

df_Q <- data.frame(doc_id = doc_ids_Q, text = text_Q, stringsAsFactors = FALSE)
corpus_Q <- Corpus(DataframeSource(df_Q))

# Convertir le texte en minuscule
corpus_Q <- tm_map(corpus_Q, content_transformer(tolower))

# Supprimer les nombres
corpus_Q <- tm_map(corpus_Q, removeNumbers)

# Supprimer les mots vides français

corpus_Q <- tm_map(corpus_Q, removeWords, stopwords("fr"))

# Supprimer les ponctuations
corpus_Q <- tm_map(corpus_Q, removePunctuation)

# Supprimer les espaces vides supplémentaires
corpus_Q <- tm_map(corpus_Q, stripWhitespace)

# Supprimer votre propre liste de mots non désirés
corpus_Q <- tm_map(corpus_Q, removeWords, c("questions","manquent","injustes","mettre",
                                            "dont","sud","semble","autres","vaut","verts",
                                            "inrae","ird","faire","autre","cas","chaque",
                                            "arrivait","mette","toutes","compte","articule",
                                            "oui","doit","individuel","voir","annuel",
                                            "entretien","ademe","démarche","pluie","espaces","sincérité",
                                            "prendre","faire","être","ceux","dom","site","incomplètes",
                                            "billets","dom","plus","envoi","donne","toutes","complexes",
                                            "cidessus","grand","avant","rendre","lieu","coût","mesures",
                                            "interessant","parait","certaines","faut","fait","etc","cirad", "Vs","tout","très","non","comme",
                                            "ça","donc","quand","aussi","bien","bonne","sait","où")) 

# Text stemming
#corpus_Q1 <- tm_map(corpus_Q1, stemDocument)
dtm_Q <- TermDocumentMatrix(corpus_Q)

#dtm_Q<- TermDocumentMatrix(corpus_Q)
dtm_Q<- weightTfIdf(dtm_Q, normalize=TRUE)
m_Q <- as.matrix(dtm_Q)
#m_Q<- t(m_Q)
v <- sort(rowSums(m_Q),decreasing=TRUE)
d <- data.frame(word = names(v),freq=v)
head(d, 200)

# nuage de mot

set.seed(1234)
wordcloud(words = d$word, freq = d$freq, min.freq = 0.5,
          max.words=200, random.order=FALSE, rot.per=0.35, 
          colors=brewer.pal(8, "Dark2"))
# Mots les plus fréquents
findFreqTerms(dtm_Q, lowfreq = 2)

barplot(d[1:20,]$freq, las = 2, names.arg = d[1:20,]$word,
        col ="lightblue", main ="Most frequent words",
        ylab = "Word frequencies")

###### text miining sur l'ensemble des questions liées à l'enquête par questionnaire ##############

# importation donne
text_E<- readLines("data/tous_commentaires.txt")
doc_ids_E <- c(1:1048)

df_E <- data.frame(doc_id = doc_ids_E, text = text_E, stringsAsFactors = FALSE)
corpus_E <- Corpus(DataframeSource(df_E))

# Convertir le texte en minuscule
corpus_E <- tm_map(corpus_E, content_transformer(tolower))

# Supprimer les nombres
corpus_E <- tm_map(corpus_E, removeNumbers)

# Supprimer les mots vides français
corpus_E <- tm_map(corpus_E, removeWords, stopwords("french"))

# Supprimer les ponctuations
corpus_E <- tm_map(corpus_E, removePunctuation)

# Supprimer les espaces vides supplémentaires
corpus_E <- tm_map(corpus_E, stripWhitespace)

# Supprimer votre propre liste de mots non désirés
corpus_E <- tm_map(corpus_E, removeWords, c("modes","mode","partie","limite","chaque","mois","bon","données","exemplarité",
                                            "devrait","pourraient","car","viens","rien","quoi","reste","jusqu",
                                            "travailler","changer","déplacement","changements","point",
                                            "avant","face","prends","effets","hors","toute","vais","effets",
                                            "toutes","bonne","choses","aussi","faut","fois","ans","tant","ceux",
                                            "autant","faudrait","autres","part","dire","sujets","recherches",
                                            "tout","peux","donc","alors","fait","fais","etc","lorsque","dont",
                                            "vie","là","vraiment","deux","notamment","où","entre","projets",
                                            "cas","ça","comme","bien","trop","plutôt","suites",
                                            "peuvent","pourrait","lié","climatiques","trouve","parce","font","chez","aime"
                                            ,"sauf","considère","réellement","veut","personnel","environnementaux","environnementales",
                                            "personnel","seule","entendu","mieux","venir","aimerais","après","commençant","pratiquant","vrai","premier"
                                            ,"hauteur","aucune","umr","pourquoi","implique","tard","retraite","certains","selon",
                                            "uniquement","semble","aller","max","chose","habite","cnrs","comment","travaux","juste",
                                            "doivent","déprimerait","ita","déplacait","peut","pieds","pense","vers","sais","inutile",
                                            "conférences","réduction","activités","sens","surtout","lien","idem","possible","tous","encore"
                                            ,"toujour","travaille","sujet","plus","tente","flics","dit","marge","teletravail","près"
                                            ,"significativement","orientée","travaillant","voir","compétents","puisque","institutsorganismesuniversités"
                                            ,"place","dès","seul","mal","aucun","limitant","soirees","dépend","doit","commun","ciradien","minime","signaux",
                                            "passer","poser","invitée","parti","privilégiant","années","essentiel","pedibus","franchi"
                                            ,"retraité","quelquefois","malheureusement","outils","également","cause","séminaires",
                                            "bcp","prend","deplacements","visioconférences","quelques","prennent","faisant","permet","forte",
                                            "veux","actuellement","voire","capitaine","faudra","mettant","compris","visvis","lesquelles","quelque",
                                            "côté","sciences","sinon","celle","septembre","certain","dépasse","mettons","sein","internationales","refusé",
                                            "veux","grâce","forte","crési","one","week","vacance","essaierai","faudra","tolteca","congres","pourquoi",
                                            "bangkok","déplacée","vue","peutêtre","chercheuse","lors","faite","etre","dois","zéro","utilisant","souhaiterais"
                                            ,"devrions","pris","individuel","falloir","base","enfin","protent","poussent",'auxquelles',"justement",
                                            "dom","-","refusent","agit","oui...","objet","objectivement","mondial","menacent","déclarations","devient",
                                            "attendrait","cirad...","allumettes","vis","infiniment","voyageant","vacances","recehrcheproportionnelle","proportionnalité",
                                            "gpl","enercoop","dernière","comùpris","réunions...","ouvrait","deviendra","call","mains","foutent","doigts",
                                            "comptent","celuici","star","posez","madrid","trabail","requete","net","lignecloud","jusquoù","deds","coute","center",
                                            "biblio","quitte","centers","cse","préférera","devraient","encontre","pis","jugeraient","httpssitesgooglecomviewtransitionscientifiquehome"
                                            ,"dis","vieux","vieille","sait","marchait","hyper","fan","expé","dessiccateurs","contrôlées","belles",
                                            "alambiqués","pdf","noel","sain","envoi","impérativement","first","xxie","soupçonne","siècle","chauffe",
                                            "nécessiterai","polémiquesgreenwashing","dan","parc","conduit","spéciaux","oiv","revues","montée","giec",
                                            "devenue","lisais","instar","uns","posées","connaissons","any","réutilise","dits","arrete","zurich","wageningen",
                                            "théorie","évidence","dîner","autour","jsutifiées","epouse","professionnelsperso","lsite",
                                            "consommée","vient","retenu","prospectionconstruction","inutilenon","fore","posé","laquelle","modère","disant",
                                            "dernièrement","new","delhi","préfèrent","interrogeais","hérésie","alibi","scientist","highlevel","jourssemaine",
                                            "constitue","importe","remplace","monte","lettre","fixé","décourageant","modulé","appelle","aboutira","exagèrent",
                                            "reconnu","posée","enceinte","plénière","acceptées","derniers","tic","connaissent","vélos","pensé","questionvoyager",
                                            "idiot","ferme","cote","contradictions","impactent","désormais","httpslabospointorgressources","meilleures","yeux","mette","représentant",
                                            "tocsin","sonnent","informent","alarment","zelande","oese","boucles","ainsi","échantillons","chacune","réalisées","qqn",
                                            "positionnait","fase","demi","ouaga","sachant","attendant","sortie","glypho","liens","process","française","devienne","simplesdéjà",
                                            "run","portes","ouvertes","franchement","divers","semblerait","marché","fixons","vienne","independent","httpswwwnaturecomarticlesz",
                                            "etait","completement","questionne","basées","enième","appelé","discussionsateliers","adaptant","zappe","taché","financières",
                                            "modifie","alleretour","emergency","problématiquesperceptions","responsabilisant","organisations","bonnemauvaise",
                                            "intitution","devais","mesri","continuent","comprennent","properso","objectifsallongement","concentrant","essais","plupart","mpl","durant",
                                            "divine","attendent","recherchedvpmt","pompagepurificationdésalinisation","intéresse","accessoirement","répondant","habillant",
                                            "aident","tech","low","focus","super","pourrions","voie","passerait","attends","tenants","répondent","quelconque","mariole",
                                            "lobbys","lobby","entrée","chauffés","résoudrait","puisqu","concourent","concourrent","actuelles","triangulaires","voulais",
                                            "souviens","disons","retrouvent","nécessitant","situ","saitpeutveut","ivoire","horssol","recherchesactivités","assurant",
                                            "arrivée","être","privées","versus","remplacent","conf","missione","voudrai","conso","certaine","ademe")) 

# Text stemming
#corpus_E <- tm_map(corpus_E, stemDocument(text_E,language = "french"))
dtm_E <- TermDocumentMatrix(corpus_E)
dtm_E<- weightTfIdf(dtm_E, normalize=TRUE)
m_E <- as.matrix(dtm_E)
v <- sort(rowSums(m_E),decreasing=TRUE)
d <- data.frame(word = names(v),freq=v)
head(d, 400)

# nuage de mot

set.seed(1234)
wordcloud(words = d$word, freq = d$freq, min.freq = 5 ,
          max.words=600, random.order=FALSE, rot.per=0.35, 
          colors=brewer.pal(8, "Dark2"))
# frequence des mots barplot

barplot(d[1:10,]$freq, las = 2, names.arg = d[1:10,]$word,
        col ="lightblue", main ="Mots peu fréquent",
        ylab = "Word frequencies")





#### Importation tableau commentaire :  analyse qualitative ##############

sheet_names <- excel_sheets(here::here("data","tableau_qualitatif.g.xlsx"))

for(name in sheet_names){
  name2 <- paste0("sheet_",name) #modification du nom de la feuille 
  data<- read_excel(here::here("data","tableau_qualitatif.g.xlsx"), sheet = name, col_names = TRUE, skip = 1)
  data$sheet <- name
  assign(name2,data)  
  #rm(name, name2,data)
}

#tab_test1<- sheet_TPC[,c(5,7,9,11,13,15,17,19,21,23,25,27,29,30)]

chi2_names<- excel_sheets(here::here("data","tableau_test_chi2.xlsx"))

for(name in chi2_names){
  name2 <- paste0("sheet_",name) #modification du nom de la feuille 
  data<- read_excel(here::here("data","tableau_test_chi2.xlsx"), sheet = name, col_names = TRUE, skip = 1)
  #data$sheet <- name
  assign(name2,data)  
  rm(name, name2,data)
}

sheet_PC<- sheet_PC[,-1]
sheet_PCD<- sheet_PCD[,-1]
sheet_PI<- sheet_PI[,-1]

rownames(sheet_PC)<- c("cluster_Pro changement","cluster_Nuances", "cluster_Réticents", "cluster_Ne se prononce pas")
rownames(sheet_PCD)<- c("cluster_Pro changement","cluster_Nuances", "cluster_Réticents", "cluster_Ne se prononce pas")
rownames(sheet_PI)<- c("cluster_Pro changement","cluster_Nuances", "cluster_Réticents", "cluster_Ne se prononce pas")

sheet_PC<- t(sheet_PC)
sheet_PCD<- t(sheet_PCD)
sheet_PI<- t(sheet_PI)

# H0 : On fait l'hypothèse qu'il y a une rélation d'indépendance entre nos clusters et les colonnes du tableau
# alpha : 5% = 0.05
test1<- chisq.test(sheet_PC) # p-value 0.03 < 0.05 significativité on rejette l'hypothèse nulle (donc il y' a une dépendance entre les clusters et les colonnes du tableau)
test2<- chisq.test(sheet_PCD) # p-value 0.12 > 0.05 on ne rejette pas l'hypothèse nulle
test3<- chisq.test(sheet_PI) # p-value 0.34 > 0.05 on ne rejette pas l'hypothèse nulle

tab<- read_excel("data/tab_synthese.xlsx", col_names = TRUE, skip = 1) %>% 
  as.data.frame()
rownames(tab)<- c("Q. Pratiques individuelles","Q. Pratiques collectives", "Q. Pratiques institutionnelles")
tab<- tab[,-1]



# importation des significations des codes

code_PC<- read_excel(here::here("data/Codes_PC.xlsx"), col_names = TRUE, skip = 1) %>% 
  as.data.frame()
