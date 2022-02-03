library("MASS")
library("jpeg")
library("dbscan")
library("Rmixmod")
library("mclust")
images<-list()
for(i in 0:4737){
  images[[i+1]]<-readJPEG(paste0("animaux_II/",i,".jpg"))
}
display(images[[1]])
dim(images[[1]][,,1])    ## 128 x 128
red<-matrix(nrow=2500,ncol=128 * 128)
for(i in 1:nrow(red))red[i,]<-c(images[[i]][,,1])
green<-matrix(nrow=2500,ncol=128 * 128)
for(i in 1:nrow(green))green[i,]<-c(images[[i]][,,2])
blue<-matrix(nrow=2500,ncol=128 * 128)
for(i in 1:nrow(blue))blue[i,]<-c(images[[i]][,,3])

mat_image<-cbind(red,green,blue)
acp_image<-prcomp(mat_image,scale. = F, center = F)


## Fonction pour visualiser une image à partir de la matrice avec toutes les images
display_image<-function(image_totale,numero){
  image <- array(image_totale[numero,], dim = c(128,128,3))
  image[image < 0] <- 0
  image[image > 1] <- 1 
  display(image)
}
k<-1000
acp_reduit<- acp_image$x[,1:k]
image_reduit <- acp_reduit %*% t(acp_image$rotation[,1:k])
## On sauvegarde les 20 premières images pour pouvoir comparer
#imk200<-image_reduit[1:20,]
#imk500<-image_reduit[1:20,]
imk1000<-image_reduit[1:20,]
par(mfrow=c(2,2))
display_image(imk200,6)
display_image(imk500,6)
display_image(imk1000,6)
display(images[[6]])
display_image(imk200,7)
display_image(imk500,7)
display_image(imk1000,7)
display(images[[7]])
display_image(imk200,8)
display_image(imk500,8)
display_image(imk1000,8)
display(images[[8]])
display_image(imk200,9)
display_image(imk500,9)
display_image(imk1000,9)
display(images[[9]])
display_image(imk200,10)
display_image(imk500,10)
display_image(imk1000,10)
display(images[[10]])
## On reconnait bien les animaux avec des points ou des stries (tigres, léopards)
## mais moins bien ceux avec un pelage uniforme comme les lions
memory.size()

dev.off()
m = 50
knn_image<-kNNdist(acp_reduit,k=m)
plot(sort(knn_image,decreasing=FALSE),type="l",xlab="Nombre de points",
     ylab="Distance au 50-ème point le plus proche",
     main="Représentation du nombre de points en abscisse ayant
     50 points plus proches que la distance en ordonnée")    

## eps 50 ou 55
eps_min <- 37
image_dbscan<-dbscan(acp_reduit,eps=eps_min,minPts=m);image_dbscan
classe1_dbscan<-which(image_dbscan$cluster==1)
classe2_dbscan<-which(image_dbscan$cluster==0)
classe3_dbscan<-which(image_dbscan$cluster==2)

## 3000 images (k=1000) 431/2569, avec tigres, loups blancs
#acp_nolion<-acp_reduit[classe2_dbscan,]
acp_classe1<-acp_reduit[classe2_dbscan,]

## Classe 1
display_image(image_reduit,classe1_dbscan[sample(seq(1,length(classe1_dbscan),by=1),1,TRUE)])

## Beaucoup de tigres oranges, léopards dans la classe 2
## Classe 2
display_image(image_reduit,classe2_dbscan[sample(seq(1,length(classe2_dbscan),by=1),1,TRUE)])

for(i in 1:length(classe1_dbscan)){
  writeJPEG(images[[classe1_dbscan[i]]],paste0("Classe1MA/",i,".jpg"))
}
for(i in 1:length(classe2_dbscan)){
  writeJPEG(images[[classe2_dbscan[i]]],paste0("Classe2MA/",i,".jpg"))
}
for(i in 1:length(classe3_dbscan)){
  writeJPEG(images[[classe3_dbscan[i]]],paste0("Classe3MA/",i,".jpg"))
}



## Grand nombre de données donc il est difficile de tester plusieurs points de départ pour les k means
kmeans_image<-kmeans(acp_reduit,centers=2,nstart=2)
table(kmeans_image$cluster)
classe1<-which(kmeans_image$cluster==1)
classe2<-which(kmeans_image$cluster==2)
classe3<-which(kmeans_image$cluster==3)
classe4<-which(kmeans_image$cluster==4)
classe5<-which(kmeans_image$cluster==5)
classe6<-which(kmeans_image$cluster==6)
classe7<-which(kmeans_image$cluster==7)
classe8<-which(kmeans_image$cluster==8)
classe9<-which(kmeans_image$cluster==9)
classe10<-which(kmeans_image$cluster==10)

for(i in 1:length(classe1)){
  writeJPEG(images[[classe1[i]]],paste0("kmeans1/",i,".jpg"))
}
for(i in 1:length(classe2)){
  writeJPEG(images[[classe2[i]]],paste0("kmeans2/",i,".jpg"))
}
for(i in 1:length(classe3)){
  writeJPEG(images[[classe3[i]]],paste0("kmeans3/",i,".jpg"))
}
for(i in 1:length(classe4)){
  writeJPEG(images[[classe4[i]]],paste0("D:/Gabriel/kmeans4/",i,".jpg"))
}
## Classe 1
display_image(image_reduit,classe1[sample(seq(1,length(classe1),by=1),1,TRUE)])

## Classe 2
display_image(image_reduit,classe2[sample(seq(1,length(classe2),by=1),1,TRUE)])

## Classe 3
display_image(image_reduit,classe3[sample(seq(1,length(classe3),by=1),1,TRUE)])

## Classe 4
display_image(image_reduit,classe4[sample(seq(1,length(classe4),by=1),1,TRUE)])

## Classe 5
display_image(image_reduit,classe5[sample(seq(1,length(classe5),by=1),1,TRUE)])

## Classe 6
display_image(image_reduit,classe6[sample(seq(1,length(classe6),by=1),1,TRUE)])

## Classe 7
display_image(image_reduit,classe7[sample(seq(1,length(classe7),by=1),1,TRUE)])

## Classe 8
display_image(image_reduit,classe8[sample(seq(1,length(classe8),by=1),1,TRUE)])

## Classe 9
display_image(image_reduit,classe9[sample(seq(1,length(classe9),by=1),1,TRUE)])

## Classe 10
display_image(image_reduit,classe10[sample(seq(1,length(classe10),by=1),1,TRUE)])



## CAH
dist_image<-dist(acp_reduit)
dist2<-dist(acp_classe1)

## Single
hclust_image_single<-hclust(dist_image,method="single")
plot(hclust_image_single) ## 2 Classes

## Ward (k=4 le 4ième concentre des animaux blancs)
hclust_image_ward<-hclust(dist_image,method="ward.D")
plot(hclust_image_ward,labels=FALSE)
clusters_ward<-rect.hclust(hclust_image_ward,k=4);clusters_ward
display_image(image_reduit,clusters_ward[[1]][sample(seq(1,length(clusters_ward[[1]]),by=1),1,TRUE)])
display_image(image_reduit,clusters_ward[[2]][sample(seq(1,length(clusters_ward[[2]]),by=1),1,TRUE)])
display_image(image_reduit,clusters_ward[[3]][sample(seq(1,length(clusters_ward[[3]]),by=1),1,TRUE)])
display_image(image_reduit,clusters_ward[[4]][sample(seq(1,length(clusters_ward[[4]]),by=1),1,TRUE)])
display_image(image_reduit,clusters_ward[[5]][sample(seq(1,length(clusters_ward[[5]]),by=1),1,TRUE)])
display_image(image_reduit,clusters_ward[[6]][sample(seq(1,length(clusters_ward[[6]]),by=1),1,TRUE)])

hclust_class1<-hclust(dist_image,method="ward.D")
plot(hclust_class1,labels=FALSE)
clusters_class1<-rect.hclust(hclust_class1,k=5)
#clusters_class1<-rect.hclust(hclust_class1,k=2)
for(i in 1:length(clusters_class1[[1]])){
  writeJPEG(images[[clusters_class1[[1]][i]]],paste0("Classe2CAH/",i,".jpg"))
}
for(i in 1:length(clusters_class1[[2]])){
  writeJPEG(images[[clusters_class1[[2]][i]]],paste0("Classe3CAH/",i,".jpg"))
}
for(i in 1:length(clusters_class1[[3]])){
  writeJPEG(images[[clusters_class1[[3]][i]]],paste0("Classe4CAH/",i,".jpg"))
}
for(i in 1:length(clusters_class1[[4]])){
  writeJPEG(images[[clusters_class1[[4]][i]]],paste0("Classe5CAH/",i,".jpg"))
}
for(i in 1:length(clusters_class1[[5]])){
  writeJPEG(images[[clusters_class1[[5]][i]]],paste0("Classe6CAH/",i,".jpg"))
}
length(clusters_class1[[1]])
length(clusters_class1[[2]])
length(clusters_class1[[3]])
length(clusters_class1[[4]])
length(clusters_class1[[5]])

## Complete
hclust_image_complete<-hclust(dist_image,method="complete")
plot(hclust_image_complete) 

hclust2_complete<-hclust(dist2,method="complete")
plot(hclust2_complete)
clusters_complete2<-rect.hclust(hclust2_complete,k=2)
for(i in 1:length(clusters_complete2[[1]])){
  writeJPEG(images[[clusters_complete2[[1]][i]]],paste0("Classe2CAH/",i,".jpg"))
}
for(i in 1:length(clusters_complete2[[2]])){
  writeJPEG(images[[clusters_complete2[[2]][i]]],paste0("Classe3CAH/",i,".jpg"))
}
for(i in 1:length(clusters_complete2[[3]])){
  writeJPEG(images[[clusters_complete2[[3]][i]]],paste0("Classe4CAH/",i,".jpg"))
}
for(i in 1:length(clusters_complete2[[4]])){
  writeJPEG(images[[clusters_complete2[[4]][i]]],paste0("Classe5CAH/",i,".jpg"))
}
## Average
hclust_image_average<-hclust(dist_image,method="average")
plot(hclust_image_average) 

mc_c2<-Mclust(acp_classe1,G=2)
summary(mc_c2)
plot(mc_c2,what="BIC")

mc_c3<-Mclust(acp_classe1,G=3)
summary(mc_c3)
plot(mc_c3,what="BIC")

mc_c4<-Mclust(acp_classe1,G=4)
summary(mc_c4)
plot(mc_c4,what="BIC")
mc1<-which(mc_c3$classification==1)
mc2<-which(mc_c3$classification==2)
mc3<-which(mc_c3$classification==3)
for(i in 1:length(mc1)){
  writeJPEG(images[[mc1[i]]],paste0("EM1/",i,".jpg"))
}
for(i in 1:length(mc2)){
  writeJPEG(images[[mc2[i]]],paste0("EM2/",i,".jpg"))
}
for(i in 1:length(mc3)){
  writeJPEG(images[[mc3[i]]],paste0("EM3/",i,".jpg"))
}


mc_c5<-Mclust(acp_classe1,G=5)
summary(mc_c5)
plot(mc_c5,what="BIC")






