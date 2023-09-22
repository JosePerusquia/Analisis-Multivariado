############################################################################
# Análisis de componentes principales     
# Autor: José A. Perusquía Cortés
############################################################################

############################################################################
# Libraries
############################################################################
library(ggplot2)
library(GGally)
library(ggthemes)
library(mvtnorm)
library(car)
library(factoextra)
library(here)
library(ggfortify)

############################################################################
# Simulamos datos
############################################################################
mu=c(0,0)
sigma=matrix(c(2,3,3,5),byrow = T,nrow=2)
set.seed(31415)
X=data.frame(rmvnorm(100,mu,sigma))


p=ggplot(X,aes(x=X1,y=X2))+
  geom_point()+
  geom_vline(xintercept=0,linewidth=.15)+
  geom_hline(yintercept=0,linewidth=.15)+
  labs(x=expression(X[1]),y=expression(X[2]))+
  theme_light()
p


############################################################################
# PCA con prcomp
############################################################################

# Con matriz de varianza y covarianza
pca=prcomp(X,scale=F)
summary(pca)

# Primer componente
x_seq=seq(-3*pca$sdev[1],3*pca$sdev[1],by=.1)
y1_seq=numeric(length(x_seq))
y2_seq=numeric(length(x_seq))

for(i in 1:length(x_seq)){
  y1_seq[i]=x_seq[i]*pca$rotation[1,1]
  y2_seq[i]=x_seq[i]*pca$rotation[2,1]
}

df_seq_x=data.frame("x1"=y1_seq,"x2"=y2_seq)

p=p+geom_line(data=df_seq_x,aes(x=x1,y=x2),col="red")
p

# Segundo componente
x_seq=seq(-3*pca$sdev[2],3*pca$sdev[2],by=.1)
y1_seq=numeric(length(x_seq))
y2_seq=numeric(length(x_seq))

for(i in 1:length(x_seq)){
  y1_seq[i]=x_seq[i]*pca$rotation[1,2]
  y2_seq[i]=x_seq[i]*pca$rotation[2,2]
}

df_seq_x=data.frame("x1"=y1_seq,"x2"=y2_seq)


p=p+geom_line(data=df_seq_x,aes(x=x1,y=x2),col="red")
p

############################################################################
# Elipse definida por los componentes principales 
############################################################################
radius=c(.5,1,1.5,2,2.5,3)

for(i in 1:6){
  ConfReg=data.frame(ellipse(mu, shape=sigma, radius=radius[i], col="red",
                             add=F,draw=F, lty=1,segments=100))
  p=p+geom_point(data=ConfReg,aes(x=x,y=y),col="blue",size=.1)
}

p

############################################################################
# Nuevas variables
############################################################################

p=ggplot(data.frame(pca$x),aes(x=PC1,y=PC2))+
  geom_point()+
  geom_vline(xintercept=0,linewidth=.25)+
  geom_hline(yintercept=0,linewidth=.25)+
  xlim(-9.5,9.5)+
  labs(x=expression(PC[1]),y=expression(PC[2]))+
  theme_light()
p

############################################################################
# Ejemplo: Calificaciones
############################################################################
calificaciones=read.table(here("../Datos/calificaciones.txt"),header = T)
str(calificaciones)

# PCA con covarianzas
pca_calificaciones=prcomp(calificaciones)
pca_calificaciones$sdev^2

# Varianza de los componentes
x=c(1,2,3,4,5)
y=pca_calificaciones$sdev^2
exp_var=data.frame("PC"=x,"Var"=y)

p=ggplot(data=exp_var,aes(x=PC,y=Var))+
  geom_col(fill="navyblue",alpha=.85)+ 
  theme_light()
p

# Biplot
biplot(pca_calificaciones, col = c( "black","steelblue"), 
       cex = c(0.8, 0.65),xlim=c(-.35,.35))


# PCA con correlaciones
pca_calificaciones=prcomp(calificaciones,scale=T)
pca_calificaciones$sdev
pca_calificaciones$rot

# Biplot
fviz_pca_biplot(pca_calificaciones,title="",
                ggtheme = theme_minimal(),geom="text")


############################################################################
# Ejemplo 2: Pizza
############################################################################
pizza=read.csv(here("../Datos/Pizza.csv"),header = T)
str(pizza)

# PCA con matriz de correlaciones
pca_pizza=prcomp(pizza[,-c(1,2)],scale=T)
summary(pca_pizza)

pcs_pizza=data.frame("PC1"=pca_pizza$x[,1],"PC2"=pca_pizza$x[,2],
                     "Brand"=pizza[,1])

# Varianza de los componentes
x=c(1,2,3,4,5,6,7)
y=pca_pizza$sdev^2
exp_var=data.frame("PC"=x,"Var"=y)

p=ggplot(data=exp_var,aes(x=PC,y=Var))+
  geom_col(fill="navyblue",alpha=.85)+ 
  theme_light()
p

# Biplot
fviz_pca_biplot(pca_pizza,title="",habillage=pizza$brand,
                ggtheme = theme_minimal(),geom="point")



