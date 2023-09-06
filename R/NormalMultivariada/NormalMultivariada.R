############################################################################
# Distribución normal multivariada y distribuciones asociadas
# Autor: José A. Perusquía Cortés
############################################################################

############################################################################
# Librerías   

library(ggplot2)
library(GGally)
library(ggthemes)
library(mvtnorm)
library(scatterplot3d)
library(expm)
library(plotly)
library(MVN)
library(car)
library(ICSNP)  # Pruebas de hipótesis Hotelling
library(here)
source(here('pruebasGraficasNormMult.R'))
############################################################################

############################################################################
# Ejemplo de distribución normal multivariada no degenerada                                               

# Parámetros
mu=c(0,0)
sigma=matrix(c(3,1,1,3),byrow = T,nrow=2)

# Muestra aleatoria y su densidad para hacer un diagrama de dispersión
multnorm.sample=rmvnorm(1000,mu,sigma)
dens=dmvnorm(multnorm.sample,mean=mu,sigma=sigma)

scatterplot3d(cbind(multnorm.sample,dens),    
              color="blue", pch="", type = "h",             
              xlab = expression(x[1]), ylab = expression(x[2]), 
              zlab = expression(f(x)))

# Densidad
x=seq(-6, 6 , length.out = 100)
y=seq(-6 ,6, length.out = 100)
grid=expand.grid(x,y)
mvds <- dmvnorm(grid,mu,sigma)
matrix_mvds <-  matrix(mvds, nrow = 100)

# Perspectiva
persp(matrix_mvds, theta =270, phi = 30, expand = 0.6, 
      shade = 0.2, col = "lightblue", xlab = expression(x),
      ylab = expression(y), zlab = expression(f))


# Transformación a distribución normal multivariada estándar
X=rmvnorm(10000,mu,sigma)
sigma_inv=solve(sigma)

ggpairs(as.data.frame(X),upper = list(continuous = wrap("cor")))

Y=matrix(nrow=10000,ncol=2)
for(i in 1:10000){
  Y[i,]=sqrtm(sigma_inv)%*%(X[i,]-mu)
}

ggpairs(as.data.frame(Y),upper = list(continuous = wrap("cor")))

# Curvas de nivel v1

z = matrix(0,nrow=100,ncol=100)

for (i in 1:100) {
  for (j in 1:100) {
    z[i,j] <- dmvnorm(c(x[i],y[j]), mean=mu,sigma=sigma)
  }
}

contour(x,y,z,xlim=c(-4,4),ylim=c(-4,4))

# Curvas de nivel v2
plot_ly(x=x,y=y,z=z,type = "contour")
############################################################################

############################################################################
# Checar normalidad. Pruebas gráficas y pruebas de Mardia, Royston y
# Henze-Zirkler

# Ejemplo datos normales
pruebas_normalidad_univariadas(multnorm.sample)
prueba_forma_cuad(multnorm.sample)

mvn(multnorm.sample,mvnTest="hz",
    multivariatePlot = "none")$multivariateNormality
mvn(multnorm.sample,mvnTest="royston",
    multivariatePlot = "none")$multivariateNormality
mvn(multnorm.sample,mvnTest="mardia",
    multivariatePlot = "none")$multivariateNormality

# Iris
pruebas_normalidad_univariadas(as.matrix(iris[,-5]))
prueba_forma_cuad(as.matrix(iris[,-5]))

mvn(as.matrix(iris[,-5]),mvnTest="hz",
    multivariatePlot = "none")$multivariateNormality
mvn(as.matrix(iris[,-5]),mvnTest="royston",
    multivariatePlot = "none")$multivariateNormality
mvn(as.matrix(iris[,-5]),mvnTest="mardia",
    multivariatePlot = "none")$multivariateNormality

###########################################################################


############################################################################
# Distribución Wishart              
set.seed(3141592)
W=rWishart(n=4,df=2,Sigma=diag(1,2))

ellipse(c(0, 0), shape=W[,,1], radius=1, col="red", lty=2,add=F,xlim=c(-2,2),
        ylim=c(-2,2))
ellipse(c(0, 0), shape=W[,,2], radius=1, col="blue", lty=2)
ellipse(c(0, 0), shape=W[,,3], radius=1, col="green", lty=2)
ellipse(c(0, 0), shape=W[,,4], radius=1, col="purple", lty=2)
############################################################################


############################################################################
# Pruebas de hipótesis                                      
mu=c(64.1,64.7)
sigma=matrix(c(191,155.6,155.6,313.5),byrow = T,nrow=2)

set.seed(32131)
multnorm.sample=data.frame(rmvnorm(203,mu,sigma))

# Sigma conocida
p=ggplot(multnorm.sample,aes(x=X1,y=X2))+
  geom_point()+
  geom_point(aes(x=60,y=60),col="red",shape=3)+
  labs(x=expression(X[1]),y=expression(X[2]))+
  theme_light()
p

# Sample mean
mu_bar=apply(multnorm.sample,2,mean);mu_bar

# Statistic
203*t(mu_bar-60)%*%solve(sigma)%*%(mu_bar-60)

# Critical value at alpha=.05
qchisq(.95,2)

# Confidence region
ConfReg=data.frame(ellipse(mu_bar, shape=sigma, 
                           radius=sqrt(qchisq(.95,2)/203), 
                           col="red",add=F,draw=F, lty=2))

p=p+geom_point(data=ConfReg,aes(x=x,y=y),col="blue",size=.1)
p

# Sigma desconocida
mu_bar=apply(multnorm.sample,2,mean)
sigma_hat=cov(multnorm.sample);sigma_hat

# Statistic
((201*203)/(2*202))*t(mu_bar-60)%*%solve(sigma_hat)%*%(mu_bar-60)

# ICSNP 
HotellingsT2(multnorm.sample,mu=c(60,60))

# Critical value at alpha=.05
qf(.99,2,201)

# Confidence region
ConfReg=data.frame(ellipse(mu_bar, shape=sigma_hat, 
                           radius=sqrt((2*202)/(201*203)*qf(.95,2,201)), 
                           col="red",add=F,draw=F, lty=2))

p+geom_point(data=ConfReg,aes(x=x,y=y),col="blue",size=.1)

############################################################################