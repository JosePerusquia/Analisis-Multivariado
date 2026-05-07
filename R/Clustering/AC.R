##################################################################
# Cluster analysis                                
# Author: Jose Antonio Perusquia Cortes
# Afil: Facultad de Ciencias - UNAM
# Module: Multivariate Analysis
##################################################################

##################################################################
# Required libraries                                            
library(cluster)            # Version 2.1.8.2
library(ggplot2)            # Version 4.0.2
library(ggthemes)           # Version 5.2.0
library(factoextra)         # Version 2.0.0
library(TeachingDemos)      # Version 2.13
library(ggdendro)           # Version 0.2.0
library(dendextend)         # Version 1.19.1
library(here)               # Version 1.0.2
##################################################################

##################################################################
# Source auxiliary functions
source(here('R/Clustering/plotsAC.R'))
##################################################################

##################################################################
# Iris data
iris_scaled = scale(iris[,-5])

# Euclidean distance matrix required for silhouette analysis
# and hierarchical clustering methods
dd = dist(iris_scaled, method = "euclidean")

# Manhattan distance matrix used for silhouette analysis
dd_man = dist(iris_scaled, method = "manhattan")
##################################################################

##################################################################
# Hierarchical clustering method (agnes)                                                   

# Agnes alternative function using hclust
hc_single = hclust(dd,method="single")

# Dendrogram
p_single = plot_dendrogram(hc_single,title='Single linkage')

# Zoom into cluster 1 region (based on dendrogram structure)
p_single + coord_cartesian(xlim=c(0,95.8),ylim=c(0,1))

# Zoom into cluster 2 region (based on dendrogram structure)
p_single + coord_cartesian(xlim=c(103,150),ylim=c(0,1.35))

# It looks like 2 clusters can be considered so we cut the
# dendrogram to have k=2
clusters_single = cutree(hc_single, k = 2)

# Silhouette analysis 
sil_single = plot_silhouette(clusters_single,dd,
                             col = c("red", "green"))

# Observations with negative silhouette widths,
# indicating possible poor cluster assignment
sil_single$MC

# Agnes with Ward method
hc_ward = hclust(dd, method = "ward.D2")

# Dendrogram
p_ward = plot_dendrogram(hc_ward,title='Ward')

# Zoom into cluster 1 region (based on dendrogram structure)
p_ward+coord_cartesian(xlim=c(0,47.35),ylim=c(0,6.5))

# Zoom into cluster 2 region (based on dendrogram structure)
p_ward+coord_cartesian(xlim=c(51,78),ylim=c(0,4.5))

# Zoom into cluster 3 region (based on dendrogram structure)
p_ward+coord_cartesian(xlim=c(82.6,147.5),ylim=c(0,8))

# It looks like 2 clusters can be considered so we cut the
# dendrogram to have k=2
clusters_ward = cutree(hc_ward, k = 3)

# Silhouette analysis 
sil_ward = plot_silhouette(clusters_ward,dd,
                      col = c("red", "green","blue"))

# Observations with negative silhouette widths,
# indicating possible poor cluster assignment
sil_ward$MC
##################################################################

##################################################################
# Divisive clustering method (diana)
hc_diana = diana(iris_scaled,metric="euclidean")

# Dendrogram
p_diana = plot_dendrogram(hc_diana,title='Diana')

# Zoom into cluster 1 region (based on dendrogram structure)
p_diana+coord_cartesian(xlim=c(0,48),ylim=c(0,5.5))

# Zoom into cluster 2 region (based on dendrogram structure)
p_diana+coord_cartesian(xlim=c(52.5,92.5),ylim=c(0,4))

# Zoom into cluster 3 region (based on dendrogram structure)
p_diana+coord_cartesian(xlim=c(97,148),ylim=c(0,3.5))

# Cut the dendrogram to have three groups
clusters_diana = cutree(as.hclust(hc_diana), k = 3)

# Silhouette analysis
sil_diana = plot_silhouette(clusters_diana,dd,
                           col = c("red", "green","blue"))

# Observations with negative silhouette widths,
# indicating possible poor cluster assignment
sil_diana$MC
##################################################################

##################################################################
# Iris and k-means

# We plot the data using PCA first with the true labels
pca_iris=prcomp(iris_scaled)
iris_true=data.frame(X=pca_iris$x[,1],Y=pca_iris$x[,2],Col=iris[,5])

p_iris=ggplot(data=iris_true,aes(x=X,y=Y,col=Col))+
       geom_point(show.legend=FALSE)+
       theme_light()+
       labs(x="PC1",y="PC2")
print(p_iris)

# K-means, plotting the results using the same PCA representation
set.seed(3141592)
iris_k_means = kmeans(iris_scaled,centers=3,iter.max=100)
iris_k_means_df = data.frame(X=pca_iris$x[,1],Y=pca_iris$x[,2],
                       Col=as.factor(iris_k_means$cluster))

p_iris_k = ggplot(data=iris_k_means_df,aes(x=X,y=Y,col=Col))+
           geom_point(show.legend=FALSE)+
           theme_light()+
           labs(x="PC1",y="PC2")
print(p_iris_k)

# We add the centroids into the plot
centroids = project_to_pca(iris_k_means$centers, pca_iris)
df = data.frame(X=centroids[,1],Y=centroids[,2])
p_iris_k+geom_point(data=df,aes(x=X,y=Y),colour="black")

# Clusters 
clusters_kmeans = iris_k_means$cluster

# Silhouette (requires Euclidean distance matrix)
sil_kmeans = plot_silhouette(clusters_kmeans,dd,
                             col = c("red", "green","blue"))

# Observations with negative silhouette widths,
# indicating possible poor cluster assignment
sil_kmeans$MC
##################################################################

##################################################################
# K-medoids also known as pam algorithm using Manhattan distance
iris_pam = pam(iris_scaled,k=3,metric="manhattan")

# Plot the groups
iris_pam_df = data.frame(X=pca_iris$x[,1],Y=pca_iris$x[,2],
                    Col=as.factor(iris_pam$clustering))
p_pam=ggplot(data=iris_pam_df,aes(x=X,y=Y,col=Col))+
  geom_point(show.legend=FALSE)+
  theme_light()+
  labs(x="PC1",y="PC2")
print(p_pam)

# Add the medoids
medoids = project_to_pca(iris_pam$medoids, pca_iris)
df = data.frame(X=medoids[,1],Y=medoids[,2])
p_pam+geom_point(data=df,aes(x=X,y=Y),colour="black")

# Clusters
clusters_pam = iris_pam$clustering

# Silhouette analysis using Manhattan distance matrix dd_man
sil_pam=plot_silhouette(clusters_pam,dd_man,
                        col = c("red", "green","blue"))

# Observations with negative silhouette widths,
# indicating possible poor cluster assignment
sil_pam$MC
##################################################################

##################################################################
# Iris multi-class membership Fanny method
iris_fanny=fanny(iris_scaled,k=3,metric="euclidean")
iris_fanny$membership

# Hard clustering
clusters_fanny = iris_fanny$clustering

# Silhouette analysis
sil_fanny = plot_silhouette(clusters_fanny,dd,
                            col = c("red", "green","blue"))

# Observations with negative silhouette widths,
# indicating possible poor cluster assignment
sil_fanny$MC
##################################################################

