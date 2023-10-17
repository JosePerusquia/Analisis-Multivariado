############################################################################
# Escalamiento multidimensional métrico
# Autor: José A. Perusquía Cortés
############################################################################

########################################################################
# Librerías
library(usmap)
library(vegan)
library(xtable)
library(MASS)
library(plot3D)
library(scatterplot3d)
library(ggplot2)
library(ggthemes)
library(geodist)
library(RColorBrewer)
library(here)
########################################################################

########################################################################
# Ejemplo1: Ciudades de Estados Unidos
p=plot_usmap(regions = "states",exclud=c("AK","HI")) + 
  theme(panel.background=element_blank())+
  geom_point(aes(x=1448560,y=-1078800),cex=.5,color="red")+
  geom_point(aes(x=985952,y=-274372),cex=.5,color="red")+
  geom_point(aes(x=-413548,y=-585116),cex=.5,color="red")+
  geom_point(aes(x=457604,y=-1683698),cex=.5,color="red")+
  geom_point(aes(x=-1640522,y=-1032437),cex=.5,color="red")+
  geom_point(aes(x=1956170,y=-1891918),cex=.5,color="red")+
  geom_point(aes(x=2156170,y=-120000),cex=.5,color="red")+
  geom_point(aes(x=-1940522,y=-532437),cex=.5,color="red")+
  geom_point(aes(x=-1640522,y=592437),cex=.5,color="red")+
  geom_point(aes(x=1956170,y=-400000),cex=.5,color="red")
p

# Matriz de distancias
a=numeric(10)
b=c(587,rep(0,9))
c=c(1212,920,rep(0,8))
d=c(701,940,879,rep(0,7))
e=c(1936,1745,831,1374,rep(0,6))
f=c(604,1188,1726,968,2339,rep(0,5))
g=c(748,713,1631,1420,2451,1092,rep(0,4))
h=c(2139,1858,949,1645,347,2594,2571,rep(0,3))
i=c(2182,1737,1021,1891,959,2734,2408,678,rep(0,2))
j=c(543,597,1494,1220,2300,923,205,2442,2329,rep(0,1))

D=rbind(a,b,c,d,e,f,g,h,i,j);D
D=D+t(D);D

# Obtenemos la solucion
n=10

# Primer paso: encontramos la matriz A y la matriz doblemente centrada B
A=-D^2/2
I=diag(x=1,n,n)
U=rep(1,n)
H=I-(U%*%t(U)/n)
B=H%*%A%*%H;B

# Segundo paso: encontrar la descomposición espectral de B
p=eigen(B)$values;p
v=eigen(B)$vectors;v


# Tercer paso: Construimos configuración en 6 dimensiones
p1=p[1:6];p1
v1=v[,(1:6)];v1

y=matrix(nrow=n,ncol=6)
for(i in 1:6){
  y[,i]=sqrt(p1[i])*v1[,i]
}

# Cuarto paso: Graficamos las priemras dos coordenadas
coords=data.frame("Lat"=y[,1],"Lon"=y[,2])
labels=c("Atl","Chic","Denv","Hous","LA","Mia","NY","SF","Seat","Wash")

p=ggplot(data=coords,aes(x=Lat,y=Lon,))+
  geom_text(show.legend=F,data=coords,
            aes(label = labels))+
  theme_light()
p

# Quinto paso (opcional): Rotamos la solución de ser requerido
coords=data.frame("Lat"=-y[,1],"Lon"=-y[,2])

p=ggplot(data=coords,aes(x=Lat,y=Lon,))+
  geom_text(show.legend=F,data=coords,
            aes(label = labels))+
  theme_light()
p

# Si no conocemos las ciudades solo graficamos los puntos
p=ggplot(data=coords,aes(x=Lat,y=Lon))+
  geom_point()+
  theme_light()
p
########################################################################

########################################################################
# Ejemplo 2: Calificaciones
calificaciones=read.table(here("../Datos/calificaciones.txt"),header = T)

# MDS y PCA
W=as.matrix(sweep(calificaciones,2,colMeans(calificaciones)))

B=W%*%t(W)
S=t(W)%*%W

Lambda_B=diag(eigen(B)$values[c(1:5)],nrow=5);Lambda_B
Lambda_S=diag(svd(S)$d,nrow=5);Lambda_S

# Coordenadas principales
Lambda_B=sqrt(Lambda_B);Lambda_B
U=eigen(B)$vec[,c(1:5)];U

pcoa=U%*%Lambda_B
head(pcoa[,c(1:2)])

# Componentes principales
pca=prcomp(calificaciones,center=T)
head(pca$x[,c(1:2)])

########################################################################

########################################################################
# Constante aditiva
mod=cmdscale(D,k=9,add=T,eig=T)
mod$eig
mod$ac

#Primeros dos componentes
coords=data.frame("Lat"=-mod$points[,1],"Lon"=-mod$points[,2])
labels=c("Atl","Chic","Denv","Hous","LA","Mia","NY","SF","Seat","Wash")

p=ggplot(data=coords,aes(x=Lat,y=Lon,))+
  geom_text(show.legend=F,data=coords,
            aes(label = labels))+
  theme_light()
p

########################################################################

########################################################################
# Non-metric MDS

# Ejemplo 3: Rollo Suizo
set.seed(314)
u=runif(1000)
v=runif(1000)

x=.5*v*sin(4*pi*v)
y=u-.5
z=.5*v*cos(4*pi*v)

X=cbind(x,y,z)

scatter3D(x,y,z,phi=180,theta=180,colkey= FALSE,col=jet.col(100))


# Escalamiento multidimiensional clásico
D=dist(X)
mod=cmdscale(D,add=F,eig=T)

#Primeros dos componentes
coords=data.frame(x=mod$points[,1],y=mod$points[,2],z=z)

ggplot(data=coords,aes(x=x,y=y))+
  geom_point(aes(colour=z))+
  scale_colour_gradientn(colours=jet.col(100))+
  theme_minimal()

# ISOMDS
iso=isoMDS(D,k=2,maxit=300)
Z1=data.frame(x=iso$points[,1],y=iso$points[,2],z=z)

ggplot(data=Z1,aes(x=x,y=y))+
  geom_point(aes(colour=z))+
  scale_colour_gradientn(colours=jet.col(100))+
  theme_minimal()

# Sammon
sam=sammon(d=D,k=2,niter=200,tol=1e-10,magic=.3)
Z3=data.frame(x=sam$points[,1],y=sam$points[,2],z=z)

ggplot(data=Z3,aes(x=x,y=y))+
  geom_point(aes(colour=z))+
  scale_colour_gradientn(colours=jet.col(100))+
  theme_minimal()

#ISOMAP
iso2=isomap(D,ndim=2,k=12)
Z2=data.frame(x=iso2$points[,1],y=iso2$points[,2],z=z)

ggplot(data=Z2,aes(x=x,y=y))+
  geom_point(aes(colour=z))+
  scale_colour_gradientn(colours=jet.col(100))+
  theme_minimal()
########################################################################
