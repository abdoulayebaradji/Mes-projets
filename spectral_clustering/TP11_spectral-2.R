# Partitionnement spectral
require(igraph)

# Exercice 1
# Simulations données
n <- 100
x <- matrix(0,nr=n,nc=2)

# On génère un vecteur de taille n contenant le numéro de la composante de chaque observation
comp <- sample(1:3,n,replace = T)
for (i in 1:n)
{
  if (comp[i]==1) # simulation pour les données de la composante 1
  {
    x[i,] <- c(rnorm(1,-3,sqrt(0.2)),rnorm(1,0,sqrt(0.1)))
  }else if (comp[i]==2) # simulation pour les données de la composante 2
  {
    x[i,] <- c(rnorm(1,0,sqrt(0.2)),rnorm(1,1,sqrt(0.1)))
  }else{ # simulation pour les données de la composante 3
    x[i,] <- c(rnorm(1,2,sqrt(0.2)),rnorm(1,-1,sqrt(0.1)))
  }
}

# Nuage de points, par la fonction plot de R ou par le package ggplot
plot(x,type="p")
require(ggplot2)
ggplot(data=data.frame(x),aes(x[,1],x[,2])) + geom_point()


# Matrice de similarité
d <- as.matrix(dist(x))
s <- exp(-d^2/(2*1^2))

# Graphe dense
p <- graph.adjacency(s,weighted=T, mode = "undirected")
p <- simplify(p)
plot(p)
E(p)$width <- E(p)$weight
plot(p)

Dcomp <- apply(s,1,sum) # diagonale de la matrice des degrés
mean(Dcomp) # degré moyen
Wcomp <- s

# Graphe du epsilon-voisinage
seuil <- quantile(s,0.65)
s_eps <- s
s_eps[s<seuil] <- 0
s_eps[s>=seuil] <- 1 # (graphe binaire)
p <- graph.adjacency(s_eps,weighted=T, mode = "undirected")
p <- simplify(p)
plot(p)
E(p)$width <- E(p)$weight
plot(p)

Deps <- apply(s_eps,1,sum) # matrice des degrés
mean(Deps) # degré moyen
Weps <- s_eps

# Graphe des kppv
k <- 2*floor(log(n))

# On construit d'abord une matrice contenant les indices des (k+1) plus proches voisins de chaque point
# On garde (k+1) voisins car le plus proche sera forcément le point lui-même
knn_indices <- matrix(0,nr=n,nc=(k+1))
for (i in 1:n)
{
  knn_indices[i,] <- order(d[i,])[1:(k+1)]
}

# Pour obtenir les voisins mutuels, on doit vérifier que si j est un voisin de i, alors i est un voisin de j
knn_mutual <- list() # on créé une liste plutôt qu'une matrice car chaque élément peut avoir une taille différent, ce qui ne sera pas bien géré par une matrice
for (i in 1:n)
{
  ind_i <- numeric() # vecteur des indices des voisins mutuels de i, c'est-à-dire des voisins de i dont i est aussi voisin
  for (j in 2:(k+1)) # on boucle sur les k ppv autres que le point lui-même, qui est dans la première colonne
  {
    ind <- knn_indices[i,j] # on récupère l'indice du j-ème voisin de i
    if (i %in% knn_indices[ind,]) # on regarde si i est dans la liste des voisins de j, c'est-à-dire s'il apparaît sur la ligne "ind" de la matrice knn_indices
    {
      ind_i <- c(ind_i,ind) # si c'est le cas, alors on a trouvé un voisin mutuel, et on ajoute cet indice au vecteur ind_i, des voisins mutuels de i
    }
  }
  knn_mutual[[i]] <- ind_i # une fois que l'on a bouclé sur tous les voisins de i, on stocke le résultat dans l'élément i de la liste knn_mutual
}

# Une fois que l'on a récupéré la liste des voisins mutuels, il suffit de récupérer les arêtes correspondantes, les autres seront mises à 0
s_knn <- 0*s # permet d'obtenir une matrice de la même taille que s, mais remplie de 0
for (i in 1:n)
{
  s_knn[i,knn_mutual[[i]]] <- s[i,knn_mutual[[i]]] # sur la ligne i de la matrice de similarité s_knn, on ne remplit que les colonnes correspondant aux voisins mutuels
}

# Graphe
p <- graph.adjacency(s_knn,weighted=T, mode = "undirected")
p <- simplify(p)
plot(p)
E(p)$width <- E(p)$weight
plot(p)

Dknn <- apply(s_knn,1,sum) # matrice des degrés
mean(Dknn) # degré moyen
Wknn <- s_knn


# Matrices Laplaciennes normalisées symétriques
Lcomp <- diag(n) - diag(1/sqrt(Dcomp)) %*% Wcomp %*%  diag(1/sqrt(Dcomp))
Leps <- diag(n) - diag(1/sqrt(Deps)) %*% Weps %*%  diag(1/sqrt(Deps))
Lknn <- diag(n) - diag(1/sqrt(Dknn)) %*% Wknn %*%  diag(1/sqrt(Dknn)) 

# Spectres des matrices : bien préciser que les matrices sont symétriques, cela permet d'optimiser les calculs
spcomp <- eigen(Lcomp,symmetric=TRUE)
speps <- eigen(Leps,symmetric=TRUE)
spknn <- eigen(Lknn,symmetric=TRUE)
# si vous obtenez une erreur ici avec spknn, il se peut que cela soit du au fait qu'un point est isolé, et que son degré est donc nul, 
# ce qui fait que la matrice n'est pas de plein rang : dans ce cas, augmentez la valeur de k

# Tracé des 10 plus petites valeurs propres. La fonction eigen renvoie les valeurs propres dans l'ordre décroissant, il faut donc récupérer les 10 dernières
ggplot(data=data.frame(values=spcomp$values[91:100],order=10:1),aes(x=order,y=values)) + geom_point() + ggtitle("Graphe complet")
ggplot(data=data.frame(values=speps$values[91:100],order=10:1),aes(x=order,y=values)) + geom_point() + ggtitle("Graphe du epsilon voisinage")
ggplot(data=data.frame(values=spknn$values[91:100],order=10:1),aes(x=order,y=values)) + geom_point() + ggtitle("Graphe des k plus proches voisins")

# Sur le graphe précédent, on identifie un saut entre la 3ème et la 4ème valeur propre, on va donc retenir 3 classes
# On calcule la matrice U des vecteurs propres associés aux 3 plus petites valeurs propres. Les vecteurs propres sont rangés dans le même ordre que
# les valeurs propres, donc il faut récupérer les 3 dernières colonnes de la matrice des vecteurs propres
Ucomp <- spcomp$vectors[,98:100]
Ueps <- speps$vectors[,98:100]
Uknn <- spknn$vectors[,98:100]

# Normalisation par lignes pour obtenir la matrice V
# Fonction de normalsiation à applisuer par ligne
rownorm <- function(v){
  return(v/sqrt(sum(v^2)))
}

# Normalisation par ligne
Vcomp <- t(apply(Ucomp,1,rownorm))
Veps <- t(apply(Ueps,1,rownorm))
Vknn <- t(apply(Uknn,1,rownorm))


# Algorithme des k-means sur les matrices obtenues précédemment
clcomp <- kmeans(Vcomp,3)
cleps <- kmeans(Veps,3)
clknn <- kmeans(Vknn,3)

# On trace les clusters identifiés par chaque méthode.
# Pour utiliser le package ggplot, on a besoin de créer un data.frame
data <- data.frame(x)
data$clcomp <- as.factor(clcomp$cluster)
data$cleps <- as.factor(cleps$cluster)
data$clknn <- as.factor(clknn$cluster)

ggplot(data=data,aes(X1,X2,color=clcomp)) + geom_point()
ggplot(data=data,aes(X1,X2,color=cleps)) + geom_point()
ggplot(data=data,aes(X1,X2,color=clknn)) + geom_point()

# Avec la fonction plot
plot(x,col=clcomp$cluster,pch=19)
plot(x,col=cleps$cluster,pch=19)
plot(x,col=clknn$cluster,pch=19)


# Exercice 2 : données IRM
require(jpeg)
img <- readJPEG("~/Documents/Enseignement/Lille/M2 ISN/AD/TP/Data/irm_small.jpeg")

# Matrice de similarité
# On transforme d'abord l'image, stockée sous forme de matrice, en vecteur
img_v <- as.vector(img)
d <- as.matrix(dist(img_v))

# On choisit sigma=0.5 mais on peut tester d'autres valeurs
W <- exp(-d^2/(2*0.5^2))

# Graphe du epsilon-voisinage
seuil <- quantile(W,0.75)
Weps <- W
Weps[W<seuil] <- 0
Weps[W>=seuil] <- 1

Deps <- apply(Weps,1,sum) # matrice des degrés
mean(Deps) # degré moyen

# Matrice Laplacienne normalisée symétrique
Leps <- diag(length(img_v)) - diag(1/sqrt(Deps)) %*% Weps %*%  diag(1/sqrt(Deps))

# Spectre
speps <- eigen(Leps,symmetric=TRUE)

# Plot
n <- length(speps$values)
ggplot(data=data.frame(values=speps$values[(n-9):n],order=10:1),aes(x=order,y=values)) + geom_point() + ggtitle("Graphe du epsilon voisinage")

# Sur le graphe précédent, on remarque un saut entre la 3ème et la 4ème valeur propre, on retient donc 3 classes
# Matrice des 3 premiers vecteurs propres
Ueps <- speps$vectors[,(n-2):n]

# Normalisation par lignes
rownorm <- function(v){
  return(v/sqrt(sum(v^2)))
}

Veps <- t(apply(Ueps,1,rownorm))

# k-means sur les lignes 
cleps <- kmeans(Veps,3)

matcl <- matrix(cleps$cluster,nr=nrow(img),byrow = F)

# Plots
image(img,col=grey((0:100)/100)) # image originale
image(matcl) # clusters

# Pour superposer les clusters les uns après les autres
cl1 <- matcl
cl1[cl1!=1] <- NA
cl2 <- matcl
cl2[cl2!=2] <- NA
cl3 <- matcl
cl3[cl3!=3] <- NA

image(img,col=grey((0:100)/100))
image(cl1,col="blue",add=T)
image(cl2,col="orange",add=T)
image(cl3,col="yellow",add=T)
