ar# TP 6 : k-means et CAH

# Exercice 1
# Génération de données

source("~/Documents/Enseignement/Lille/M2 ISN/Analyse de données/TP/datasets2.R")

d1 <- circles(200)
d2 <- difftaille(200)
d3 <- carres(200)
d4 <- losanges(200)

require(ggplot2)
d <- rbind(d1,d2,d3,d4)
d$data <- c(rep(c("Cercles","Tailles différentes","Carrés","Losanges"),each=400))

ggplot(data=d,aes(x1,x2)) + geom_point() + facet_wrap(~data,scales = "free") + theme(legend.position = "none")


# k-means
iw <- rep(0,10)
for (k in 1:10){
  clust_d1 <- kmeans(d1,centers=k, nstart = 1)
  iw[k] <- clust_d1$tot.withinss
}
plot(1:10,iw,type="b")

clust_d1 <- kmeans(d1,centers=5, nstart = 1)
res <- cbind(d1,clust=clust_d1$cluster)
ggplot(data = res,aes(x1,x2,color=as.factor(clust))) + geom_point()


d1_centres <- as.data.frame(clust_d1$centers)
ggplot(data=d1,aes(x=x1,y=x2)) + geom_point() + geom_point(data=d1_centres,col="red",size=3)


# CAH
dist_d1 <- dist(d1)
hclust_1 <- hclust(dist_d1,method="complete") 
plot(rev(hclust_1$height),type="b",xlim=c(0,20))
# critère du saut minimum "single linkage"

plot(hclust_1)
clusters <- rect.hclust(hclust_1, k=2)

cah_clust <- rep(1,nrow(d1))
cah_clust[clusters[[2]]] <- 2
d1$cah_clust <- cah_clust
ggplot(data = d1,aes(x1,x2,color=as.factor(cah_clust))) + geom_point()



# d2
iw <- rep(0,10)
for (k in 1:10){
  clust_d2 <- kmeans(d2,centers=k, nstart = 1)
  iw[k] <- clust_d2$tot.withinss
}
plot(1:10,iw,type="b")

clust_d2 <- kmeans(d2,centers=2, nstart = 1)
res <- cbind(d2,clust=clust_d2$cluster)
ggplot(data = res,aes(x1,x2,color=as.factor(clust))) + geom_point()

clust_d2 <- kmeans(d2,centers=3, nstart = 1)
res <- cbind(d2,clust=clust_d2$cluster)
ggplot(data = res,aes(x1,x2,color=as.factor(clust))) + geom_point()



dist_d2 <- dist(d2)
hclust_2 <- hclust(dist_d2,method="ward.D2")
plot(hclust_2)
clusters <- rect.hclust(hclust_2, k=2)
cah_clust <- rep(1,nrow(d2))
cah_clust[clusters[[2]]] <- 2
d2$cah_clust <- cah_clust
ggplot(data = d2,aes(x1,x2,color=as.factor(cah_clust))) + geom_point()
plot(rev(hclust_2$height),type="b",xlim=c(0,20))

# d3
clust_d3 <- kmeans(d3,centers=2, nstart = 1)
res <- cbind(d3,clust=clust_d3$cluster)
ggplot(data = res,aes(x1,x2,color=as.factor(clust))) + geom_point()

dist_d3 <- dist(d3)
hclust_3 <- hclust(dist_d3,method="average")
clusters <- rect.hclust(hclust_3, k=2)
cah_clust <- rep(1,nrow(d3))
cah_clust[clusters[[2]]] <- 2
d3$cah_clust <- cah_clust
ggplot(data = d3,aes(x1,x2,color=as.factor(cah_clust))) + geom_point()


# d4
clust_d4 <- kmeans(d4,centers=2, nstart = 1)
res <- cbind(d4,clust=clust_d4$cluster)
ggplot(data = res,aes(x1,x2,color=as.factor(clust))) + geom_point()

dist_d4 <- dist(d4)
hclust_4 <- hclust(dist_d4,method="single")
clusters <- rect.hclust(hclust_4, k=2)
cah_clust <- rep(1,nrow(d4))
cah_clust[clusters[[2]]] <- 2
d4$cah_clust <- cah_clust
ggplot(data = d4,aes(x1,x2,color=as.factor(cah_clust))) + geom_point()


# Exercice 2
tempe <- read.csv("~/Documents/Enseignement/Lille/M2 ISN/Analyse de données/TP/Data/temperatures.csv")

summary(tempe)

# k-means
ix <- rep(0,10)
for (k in 1:10){
  tempe_km <- kmeans(tempe[,-c(1,18)],centers=k,nstart=5)
  iw[k] <- tempe_km$tot.withinss
}
plot(1:10,iw,type="b")

cl2 <- kmeans(tempe[,-c(1,18)],centers=2,nstart=5)
cl3 <- kmeans(tempe[,-c(1,18)],centers=3,nstart=5)

tempe$kmeanscl2 <- cl2$cluster
tempe$kmeanscl3 <- cl3$cluster

# 2 clusters
tempe$X[tempe$kmeanscl2 == 1]
tempe$X[tempe$kmeanscl2 == 2]

# 3 clusters
tempe$X[tempe$kmeanscl3 == 1]
tempe$X[tempe$kmeanscl3 == 2]
tempe$X[tempe$kmeanscl3 == 3]



# CAH
dist_tempe <- dist(tempe[,-c(1,18)])
tempe_cah <- hclust(dist_tempe,method="ward.D2")
plot(tempe_cah,labels=tempe$X)
plot(rev(tempe_cah$height),type="b")

plot(tempe_cah,labels=tempe$X)
rect.hclust(tempe_cah,k=3)



d <- read.table("~/Documents/Enseignement/Lille/M2 ISN/Analyse de données/TP/Data/breast-cancer.txt",header=T)
head(d)

dist_d <- dist(d[,-d$Class])
cancer_cah <- hclust(dist_d,method = "complete")
plot(cancer_cah)
plot(rev(cancer_cah$height),type="b",xlim=c(0,100))

plot(cancer_cah)
cancer_2cl <- rect.hclust(cancer_cah,k=2)
d$cluster <- 1
d$cluster[cancer_2cl[[2]]] <- 2

summary(d[d$cluster==1,])
summary(d[d$cluster==2,])

d_long <- d %>% pivot_longer(cols=!cluster)
ggplot(data=d_long,aes(x=value,fill=as.factor(cluster))) + geom_bar() + facet_wrap(~name,scales="free")#facet_grid(cluster~name)
