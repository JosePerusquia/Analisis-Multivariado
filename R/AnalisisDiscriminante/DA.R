#########################################################################
# Análisis discriminante
# Autor: José A. Perusquía Cortés
#########################################################################

#########################################################################
# Librerías  
library(here)
library(mvtnorm)
library(ggplot2)
library(ggthemes)
library(MASS)
library(klaR)
#########################################################################

#########################################################################
# Ejemplo 1: Distribución normal misma matriz de covarianza

# Parámetros
mu1=c(7,3)
mu2=c(3,5)
Sigma=matrix(c(1,-1,-1,1.25),byrow = T,nrow=2)

# Observaciones
set.seed(314)
x1=rmvnorm(n=20,mean = mu1,sigma = Sigma)
x2=rmvnorm(n=20,mean = mu2,sigma = Sigma)
group=as.factor(c(rep(1,20),rep(2,20)))
df=data.frame(X=c(x1[,1],x2[,1]),Y=c(x1[,2],x2[,2]),G=group)

# Graficamos observaciones
p=ggplot(df,aes(x=X,y=Y,col=G))+
  geom_point(show.legend = F)+
  labs(x="",y="")+
  theme_light()
p

# Función discriminante lineal
Dx=function(x,pi1,pi2,mu_1,mu_2,sigma){
  lambda=solve(sigma)%*%(mu_1-mu_2)
  return((log(pi2/pi1)-lambda[1]*x+.5*lambda[1]*(mu_1[1]+mu_2[1])+
            .5*lambda[2]*(mu_1[2]+mu_2[2]))/lambda[2])  
}

x=seq(0,9,by=.1)
y=Dx(x,20,20,mu1,mu2,Sigma)
df1=data.frame(X=x,Y=y)

# Graficamos la región
p+geom_line(data=df1,aes(x=X,y=Y),col="black")
#########################################################################

#########################################################################
# Ejemplo 2: Distribución normal con matrices de covarianza diferentes

# Parámetros
mu1=c(7,3)
mu2=c(3,5)
Sigma1=matrix(c(1,-1,-1,1.25),byrow = T,nrow=2)
Sigma2=matrix(c(2,-1.25,-1.25,1),byrow = T,nrow=2)

# Observaciones
set.seed(314)
x1=rmvnorm(n=20,mean = mu1,sigma = Sigma1)
x2=rmvnorm(n=20,mean = mu2,sigma = Sigma2)
group=as.factor(c(rep(1,20),rep(2,20)))
df=data.frame(X=c(x1[,1],x2[,1]),Y=c(x1[,2],x2[,2]),G=group)

# Graficamos las observaciones
p=ggplot(df,aes(x=X,y=Y,col=G))+
  geom_point(show.legend = F)+
  labs(x="",y="")+
  theme_light()
p


# Agregamos la región 
mod=qda(df[,c(1:2)],df$G)
da_ggplot(mod,df,df$G)
#########################################################################

#########################################################################
# Ejemplo3: iris                                                      
Iris=iris[,c(2,1,5)]

mu1=apply(Iris[c(1:50),c(1,2)],2,mean);mu1
mu2=apply(Iris[c(51:100),c(1,2)],2,mean);mu2
mu3=apply(Iris[c(101:150),c(1,2)],2,mean);mu3

Sigma1=var(Iris[c(1:50),c(1,2)]);Sigma1
Sigma2=var(Iris[c(51:100),c(1,2)]);Sigma2
Sigma3=var(Iris[c(101:150),c(1,2)]);Sigma3

Spool=(49/147)*(Sigma1+Sigma2+Sigma3);Spool

p=ggplot(Iris,aes(x=Sepal.Width,y=Sepal.Length,col=Species))+
  geom_point(show.legend=F)+
  labs(x="Sepal Width",y="Sepal Length")+
  theme_minimal()
p

# Setosa vs Virginica
x=seq(2,4.5,by=.1)
y=Dx(x,20,20,mu1,mu2,Spool)
df1=data.frame(X=x,Y=y)
p+geom_line(data=df1,aes(x=X,y=Y),col="black")

# Setosa vs Versicolor
x=seq(2,4.5,by=.1)
y=Dx(x,20,20,mu1,mu3,Spool)
df1=data.frame(X=x,Y=y)
p+geom_line(data=df1,aes(x=X,y=Y),col="black")

# Versicolor vs Virginica
x=seq(2,4.5,by=.1)
y=Dx(x,20,20,mu2,mu3,Spool)
df1=data.frame(X=x,Y=y)
p+geom_line(data=df1,aes(x=X,y=Y),col="black")

# Graficamos las tres regiones simplificadas
x=seq(2,4.5,by=.1)
y=Dx(x,20,20,mu1,mu2,Spool)
df1=data.frame(X=x,Y=y)
p=p+geom_segment(data=df1,aes(x=X[1],y=Y[1],yend=6.1655,xend=3.5476),
                 col="black",linewidth=.1)

x=seq(2,4.5,by=.1)
y=Dx(x,20,20,mu2,mu3,Spool)
df2=data.frame(X=x,Y=y)
p=p+geom_segment(data=df2,aes(x=X[1],y=Y[1],yend=6.1655,xend=3.5476),
                 col="black",linewidth=.1)

x=seq(2,4.5,by=.1)
y=Dx(x,20,20,mu1,mu3,Spool)
df3=data.frame(X=x,Y=y)
p=p+geom_segment(data=df3,aes(xend=X[26],yend=Y[26],y=6.1655,x=3.5476),
                 col="black",linewidth=.1)
p
#########################################################################