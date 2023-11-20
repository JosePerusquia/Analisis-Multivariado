################################################################################
# Análisis de conglomerados
# Autor: José A. Perusquía Cortés
################################################################################

################################################################################
# Librerías                                                  
################################################################################
library(cluster)
library(ggplot2)
library(ggthemes)
library(factoextra)
library(TeachingDemos)
library(ggdendro)
library(dendextend)

################################################################################
# Ejemplo 1                                                 

# Matriz de disimilitueds
D=matrix(c(0,7,1,9,8,7,0,6,3,5,
           1,6,0,8,7,9,3,8,0,4,
           8,5,7,4,0),nrow=5,ncol=5,byrow=T);D

# Agnes con el vecino más cercano
singleLin=agnes(D,diss=T,method="single")
pltree(singleLin, cex = 1, hang = -1,main="Dendograma",sub="",xlab="")

# Agnes con el vecino más lejano produce el mismo resultado 
completeLin=agnes(D,diss=T,method="complete")
pltree(completeLin, cex = 1, hang = -1,main="Dendograma",sub="",xlab="")
################################################################################

################################################################################
# Iris con Agnes                                                     

# Creamos la matriz de distancias euclidianas
dd = dist(iris[,-5], method = "euclidean")

# Función alternativa para crear el clustering de agnes con método single
hc_s = hclust(dd,method="single")

# Creamos el dendograma
dend = as.dendrogram(hc_s)
dend_data = dendro_data(dend, type = "rectangle")

# Graficamos el dendograma completo
p = ggplot(dend_data$segments) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+
  geom_text(data = dend_data$labels, aes(x, y, label = label),
            hjust = 1, angle = 90, size = 2)+
  ylim(-.15, max(hc_s$height))+
  theme_minimal()+
  labs(x="",y="Height")
p

# Hacemos zoom al primer grupo con la opción coord_cartesian
p = ggplot(dend_data$segments) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+
  geom_text(data = dend_data$labels, aes(x, y, label = label),
            hjust = 1, angle = 90, size = 3)+
  ylim(-.15, max(hc_s$height))+
  theme_minimal()+
  labs(x="",y="Height")+
  coord_cartesian(xlim=c(0,48),ylim=c(0,.75))
p

# Segundo grupo
p = ggplot(dend_data$segments) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+
  geom_text(data = dend_data$labels, aes(x, y, label = label),
            hjust = 1, angle = 90, size = 2)+
  ylim(-.15, max(hc_s$height))+
  theme_minimal()+
  labs(x="",y="Height")+
  coord_cartesian(xlim=c(55,147.5),ylim=c(0,.8))
p

# Silhouette para analizar el clustering 
si= silhouette(cutree(hc_s, k =2), dd)
plot(si,col = c("red", "green"),main="Silhouette")
si[which(si[,3]<0),]

# Función alternativa para crear el clustering de agnes con método Ward
hc_w = hclust(dd, method = "ward.D2")

# Creamos el dendograma
dend = as.dendrogram(hc_w)
dend_data = dendro_data(dend, type = "rectangle")

# Graficamos el dendograma completo
p = ggplot(dend_data$segments) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+
  geom_text(data = dend_data$labels, aes(x, y, label = label),
            hjust = 1, angle = 90, size = 2)+
  ylim(-.15, max(hc_w$height))+
  theme_minimal()+
  labs(x="",y="Height")
p

# Hacemos zoom al primer grupo con la opción coord_cartesian
p = ggplot(dend_data$segments) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+
  geom_text(data = dend_data$labels, aes(x, y, label = label),
            hjust = 1, angle = 90, size = 3)+
  ylim(-.15, max(hc_w$height))+
  theme_minimal()+
  labs(x="",y="Height")+
  coord_cartesian(xlim=c(0,48),ylim=c(0,4.25))
p

# Segundo grupo
p = ggplot(dend_data$segments) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+
  geom_text(data = dend_data$labels, aes(x, y, label = label),
            hjust = 1, angle = 90, size = 3)+
  ylim(-.15, max(hc_w$height))+
  theme_minimal()+
  labs(x="",y="Height")+
  coord_cartesian(xlim=c(52,85),ylim=c(0,5.1))
p

# Tercer grupo
p = ggplot(dend_data$segments) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+
  geom_text(data = dend_data$labels, aes(x, y, label = label),
            hjust = 1, angle = 90, size = 3)+
  ylim(-.15, max(hc_w$height))+
  theme_minimal()+
  labs(x="",y="Height")+
  coord_cartesian(xlim=c(89.3,147.5),ylim=c(0,7))
p

# Silhouette para analizar el clustering 
si= silhouette(cutree(hc_w, k =3), dd)
plot(si,col = c("red", "green","blue"),main="Silhouette")
si[which(si[,3]<0),]
################################################################################


################################################################################
# Iris con Diana
div=diana(iris[,-5],metric="euclidean")

dend = as.dendrogram(div)
dend_data = dendro_data(dend, type = "rectangle")

# Dendograma complete
p = ggplot(dend_data$segments) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+
  geom_text(data = dend_data$labels, aes(x, y, label = label),
            hjust = 1, angle = 90, size = 2)+
  ylim(-.15, max(div$height))+
  theme_minimal()+
  labs(x="",y="Height")

p


# Primer grupo
p = ggplot(dend_data$segments) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+
  geom_text(data = dend_data$labels, aes(x, y, label = label),
            hjust = 1, angle = 90, size = 3)+
  ylim(-.15, max(div$height))+
  theme_minimal()+
  labs(x="",y="Height")+
  coord_cartesian(xlim=c(0,51),ylim=c(0,3))
p

# Segundo grupo
p = ggplot(dend_data$segments) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+
  geom_text(data = dend_data$labels, aes(x, y, label = label),
            hjust = 1, angle = 90, size = 3)+
  ylim(-.15, max(div$height))+
  theme_minimal()+
  labs(x="",y="Height")+
  coord_cartesian(xlim=c(56.25,110.75),ylim=c(0,3))

p

# Tercer grupo
p = ggplot(dend_data$segments) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+
  geom_text(data = dend_data$labels, aes(x, y, label = label),
            hjust = 1, angle = 90, size = 3)+
  ylim(-.15, max(div$height))+
  theme_minimal()+
  labs(x="",y="Height")+
  coord_cartesian(xlim=c(115,150),ylim=c(0,3))
p

## Cortamos en tres grupos
div_groups <- cutree(as.hclust(div), k = 3)

si_diana=silhouette(div_groups,dist=dd)
plot(si_diana, col = c("red", "green", "blue"),main="Silhouette")
si_diana[which(si_diana[,3]<0),]
################################################################################

############# Alternativo
  dend <- iris[,-5]  %>% dist %>% 
  hclust %>% as.dendrogram %>%
  set("branches_k_color", k=3) %>% set("branches_lwd", 1.2) %>%
  set("labels_colors") %>% set("labels_cex", c(.5)) %>% 
  set("leaves_pch", 19) %>% set("leaves_col", c("blue")) %>%
  set("leaves_cex", .5)

plot(dend)

################################################################################
# Iris y k-means

res=kmeans(iris[,-5],centers=3,iter.max=100)

pca_iris=prcomp(iris[,-5])
iris_true=data.frame(X=pca_iris$x[,1],Y=pca_iris$x[,2],
                     Col=iris[,5])

p=ggplot(data=iris_true,aes(x=X,y=Y,col=Col))+
  geom_point(show.legend=F)+
  theme_light()+
  labs(x="PC1",y="PC2")
p

iris_kmeans=data.frame(X=pca_iris$x[,1],Y=pca_iris$x[,2],
                       Col=as.factor(res$cluster))

p=ggplot(data=iris_kmeans,aes(x=X,y=Y,col=Col))+
  geom_point(show.legend=F)+
  theme_light()+
  labs(x="PC1",y="PC2")
p

centroids=matrix(nrow=3,ncol=2)
centroids[1,]=((res$centers[1,]-pca_iris$center)%*%pca_iris$rotation)[c(1:2)]
centroids[2,]=((res$centers[2,]-pca_iris$center)%*%pca_iris$rotation)[c(1:2)]
centroids[3,]=((res$centers[3,]-pca_iris$center)%*%pca_iris$rotation)[c(1:2)]
centroids

df=data.frame(X=centroids[,1],Y=centroids[,2])
p+geom_point(data=df,aes(x=X,y=Y),colour="black")

ind1=which(res$cluster==2)
ind2=which(res$cluster==3)
ind3=which(res$cluster==1)

clusters=numeric(150)
clusters[ind1]=1
clusters[ind2]=2
clusters[ind3]=3

dd = dist(iris[,-5], method = "euclidean")

si_kmeans=silhouette(res$cluster,dd)
plot(si_kmeans,col = c("red", "green", "blue"),main="Silhouette")

################################################################################

################################################################################
# medoide

res=pam(iris[,-5],k=3,metric="manhattan")

ind1=which(res$clustering==1)
ind2=which(res$clustering==3)
ind3=which(res$clustering==2)

clusters=numeric(150)
clusters[ind1]=1
clusters[ind2]=2
clusters[ind3]=3

iris_pam=data.frame(X=pca_iris$x[,1],Y=pca_iris$x[,2],
                    Col=as.factor(res$clustering))

p=ggplot(data=iris_pam,aes(x=X,y=Y,col=Col))+
  geom_point(show.legend=F)+
  theme_light()+
  labs(x="PC1",y="PC2")
p

medoids=matrix(nrow=3,ncol=2)
medoids[1,]=((res$medoids[1,]-pca_iris$center)%*%pca_iris$rotation)[c(1:2)]
medoids[2,]=((res$medoids[2,]-pca_iris$center)%*%pca_iris$rotation)[c(1:2)]
medoids[3,]=((res$medoids[3,]-pca_iris$center)%*%pca_iris$rotation)[c(1:2)]
medoids

df=data.frame(X=medoids[,1],Y=medoids[,2])
p+geom_point(data=df,aes(x=X,y=Y),colour="black")

si_pam=silhouette(res,dd)
plot(si_pam,col = c("red", "green", "blue"),main="Silhouette")
si_pam[which(si_pam[,3]<0),]

################################################################################

################################################################################
# Iris fanny
res_fanny=fanny(iris[,-5],k=3,metric="euclidean")
res_fanny$membership
res_fanny$clustering

si_fanny=silhouette(res_fanny$clustering,dd)

plot(si_fanny,col = c("red", "green", "blue"),main="Silhouette")
si_fanny[which(si_fanny[,3]<0),]
################################################################################