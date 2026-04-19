###########################################################
# Multidimensional scaling (MDS)                                         
# Author: Jose Antonio Perusquia Cortes
# Afil: Facultad de Ciencias - UNAM
# Module: Multivariate analysis
###########################################################

###########################################################
# Required libraries
library(usmap)                # Version 0.7.1
library(vegan)                # Version 2.6-10
library(xtable)               # Version 1.8-4
library(MASS)                 # Version 7.3-60
library(plot3D)               # Version 1.4.1
library(scatterplot3d)        # Version 0.3-44
library(ggplot2)              # Version 3.5.2
library(ggthemes)             # Version 5.1.0
library(geodist)              # Version 0.1.0
library(RColorBrewer)         # Version 1.1-3
library(sf)                   # Version 1.0-19
library(here)                 # Version 1.0.1
###########################################################

###########################################################
# Example: Cities of the United States
cities = read.csv(here("R/DimensionalityReduction/Data/dist_US.txt"),
                    header=F,sep=',')
# Number of cities
n = length(cities)

# Map of the US
p=plot_usmap(regions = "states",exclud=c("AK","HI")) + 
  theme(panel.background=element_blank())+
  geom_point(x=1448560,y=-1078800,cex=.5,color="red")+
  geom_point(x=985952,y=-274372,cex=.5,color="red")+
  geom_point(x=-413548,y=-585116,cex=.5,color="red")+
  geom_point(x=457604,y=-1683698,cex=.5,color="red")+
  geom_point(x=-1640522,y=-1032437,cex=.5,color="red")+
  geom_point(x=1956170,y=-1891918,cex=.5,color="red")+
  geom_point(x=2156170,y=-120000,cex=.5,color="red")+
  geom_point(x=-1940522,y=-532437,cex=.5,color="red")+
  geom_point(x=-1640522,y=592437,cex=.5,color="red")+
  geom_point(x=1956170,y=-400000,cex=.5,color="red")
p

# Distance matrix
D_cities = as.matrix(cities)

# First step: Obtain double-centred matrix B
A = -D_cities^2/2
I = diag(x=1,n,n)
U = rep(1,n)
H = I-(U%*%t(U)/n)
B = H%*%A%*%H;B

# Second step: find spectral decomposition of B, we keep
# dominant eigenvalues
eig = eigen(B)
l = eig$values
v = eig$vectors
k = length(which(l>0));k

# Third step: Since there are 6 positive eigenvalues
# we build a solution on six dimensions 
l1=l[1:k];l1
v1=v[,(1:k)];v1

y=matrix(nrow=n,ncol=k)
for(i in 1:k){
  y[,i]=sqrt(l1[i])*v1[,i]
}

# Fourth step: Plot the solution in 2 dimensions
coords=data.frame("Lat"=y[,1],"Lon"=y[,2])
labels=c("Atl","Chic","Denv","Hous","LA","Mia",
         "NY","SF","Seat","Wash")
p=ggplot(data=coords,aes(x=Lat,y=Lon,))+
  geom_text(show.legend=F,data=coords,
            aes(label = labels))+
  theme_minimal();p

# Fifth step: Rotate the solution since MDS solutions are 
# invariant to rotation, reflection, and translation
coords=data.frame("Lat"=-y[,1],"Lon"=-y[,2])
p=ggplot(data=coords,aes(x=Lat,y=Lon,))+
  geom_text(show.legend=F,data=coords,
            aes(label = labels))+
  theme_minimal();p

# In case the names of the cities are not known we just
# plot the dots
p=ggplot(data=coords,aes(x=Lat,y=Lon))+
  geom_point()+
  theme_minimal();p

# Additive constant using cmdscale function
mod = cmdscale(D_cities,add=T,eig=T)
mod$eig
mod$ac

# First two components
coords=data.frame("Lat"=-mod$points[,1],
                  "Lon"=-mod$points[,2])
labels=c("Atl","Chic","Denv","Hous","LA","Mia",
         "NY","SF","Seat","Wash")

p=ggplot(data=coords,aes(x=Lat,y=Lon,))+
  geom_text(show.legend=F,data=coords,
            aes(label = labels))+
  theme_minimal()
p
###########################################################

###########################################################
# Example: Capitals of Europe
capitals = read.csv(here("R/DimensionalityReduction/Data/dist_europa.csv"),
                    header=T,sep=',',row.names=1)
n = length(capitals)

# Distance matrix
D_capitals=as.matrix(capitals)
diag(D_capitals) = 0

# Multidimensional scaling using cmdscale function
mod=cmdscale(D_capitals,k=2,add=T,eig=T)

# Two dimensional coordinates
coords=data.frame("Lat"=mod$points[,1],
                  "Lon"=mod$points[,2])

p=ggplot(data=coords,aes(x=Lat,y=Lon,))+
  geom_point(show.legend=F,data=coords)+
  theme_light()
p
###########################################################

###########################################################
# Example: Grades
grades=read.table(here("R/DimensionalityReduction/Data/calificaciones.txt"),
                  header = T)

# Comparison between MDS and PCA, since we have a
# Euclidean matrix both methods yield the same result
# up to a rotation
W=as.matrix(sweep(grades,2,colMeans(grades)))
B=W%*%t(W)
S=t(W)%*%W

Lambda_B=diag(eigen(B)$values[c(1:5)],nrow=5);Lambda_B
Lambda_S=diag(svd(S)$d,nrow=5);Lambda_S

# Principal coordinates
Lambda_B=sqrt(Lambda_B);Lambda_B
U=eigen(B)$vec[,c(1:5)];U
pcoa=U%*%Lambda_B
head(pcoa[,c(1:2)])

# PCA 
pca=prcomp(grades,center=T)
head(pca$x[,c(1:2)])
###########################################################

###########################################################
# Non-metric MDS

# Example: Swiss roll
set.seed(314)
u=runif(1000)
v=runif(1000)

x=.5*v*sin(4*pi*v)
y=u-.5
z=.5*v*cos(4*pi*v)

X=cbind(x,y,z)
scatter3D(x,y,z,phi=180,theta=180,colkey= FALSE,
          col=jet.col(100))

# Classic mds
D=dist(X)
mod=cmdscale(D,add=F,eig=T)

# First two components
coords=data.frame(x=mod$points[,1],y=mod$points[,2],z=z)
ggplot(data=coords,aes(x=x,y=y))+
  geom_point(aes(colour=z))+
  scale_colour_gradientn(colours=jet.col(100))+
  theme_minimal()

# ISOMDS solution
iso=isoMDS(D,k=2,maxit=300)
iso$stress

Z1=data.frame(x=iso$points[,1],y=iso$points[,2],z=z)
ggplot(data=Z1,aes(x=x,y=y))+
  geom_point(aes(colour=z))+
  scale_colour_gradientn(colours=jet.col(100))+
  theme_minimal()

# Sammon mapping
sam=sammon(d=D,k=2,niter=200,tol=1e-10,magic=.3)
sam$stress

Z3=data.frame(x=sam$points[,1],y=sam$points[,2],z=z)
ggplot(data=Z3,aes(x=x,y=y))+
  geom_point(aes(colour=z))+
  scale_colour_gradientn(colours=jet.col(100))+
  theme_minimal()

# ISOMAP replaces Euclidean distances with graph-based 
# geodesic distances before applying classical MDS
iso2=isomap(D,ndim=2,k=8)
Z2=data.frame(x=iso2$points[,1],y=iso2$points[,2],z=z)

ggplot(data=Z2,aes(x=x,y=y))+
  geom_point(aes(colour=z))+
  scale_colour_gradientn(colours=jet.col(100))+
  theme_minimal()
###########################################################