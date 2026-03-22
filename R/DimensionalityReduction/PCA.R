###########################################################
# Principal components analysis (PCA)                                            
# Author: Jose Antonio Perusquia Cortes
# Afil:  Facultad de Ciencias - UNAM
# Module: Multivariate Analysis
###########################################################

###########################################################
# Required libraries 
library(ggplot2)        # Version 4.0.2
library(ggthemes)       # Version 5.2.0
library(mvtnorm)        # Version 1.3-3
library(car)            # Version 3.1-5
library(factoextra)     # Version 2.0.0
library(here)           # Version 1.0.2
library(ggfortify)      # Version 0.4.19
###########################################################

###########################################################
# Simulate data a from bivariate normal distribution  
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
print(p)

# Principal component analysis using covariance matrix
# note that by default it centers the data
pca=prcomp(X,scale.=F,center=T)
summary(pca)

# Plot first component
x_seq=seq(-3*pca$sdev[1],3*pca$sdev[1],by=.1)

center = colMeans(X)
y1_seq = center[1] + x_seq*pca$rotation[1,1]
y2_seq = center[2] + x_seq*pca$rotation[2,1]


df_seq_x=data.frame("x1"=y1_seq,"x2"=y2_seq)
p=p+geom_line(data=df_seq_x,aes(x=x1,y=x2),col="red")
print(p)

# Plot the second component
x_seq=seq(-3*pca$sdev[2],3*pca$sdev[2],by=.1)
y1_seq = center[1] + x_seq*pca$rotation[1,2]
y2_seq = center[2] + x_seq*pca$rotation[2,2]

df_seq_x=data.frame("x1"=y1_seq,"x2"=y2_seq)

p=p+geom_line(data=df_seq_x,aes(x=x1,y=x2),col="red")
print(p)

# The components represent the major axis of ellipses
radius=c(.5,1,1.5,2,2.5,3)
for(i in 1:6){
  ConfReg=data.frame(ellipse(mu, shape=sigma, 
                             radius=radius[i], col="red",
                             add=F,draw=F, lty=1,
                             segments=100))
  p=p+geom_point(data=ConfReg,aes(x=x,y=y),col="blue",
                 size=.1)
}
print(p)

# Plot the new variables
ggplot(data.frame(pca$x),aes(x=PC1,y=PC2))+
  geom_point()+
  geom_vline(xintercept=0,linewidth=.25)+
  geom_hline(yintercept=0,linewidth=.25)+
  labs(x=expression(PC[1]),y=expression(PC[2]))+
  theme_light()
###########################################################

###########################################################
# Example: Grades
grades=read.table(here("R/DimensionalityReduction/Data/calificaciones.txt"),
                  header = T)

# PCA with covariance matrix
pca_grades=prcomp(grades)
summary(pca_grades)

# Variance of the components and barplot 
x=c(1,2,3,4,5)
y=pca_grades$sdev^2;y
exp_var=data.frame("PC"=x,"Var"=y)

ggplot(data=exp_var,aes(x=PC,y=Var))+
  geom_col(fill="navyblue",alpha=.85)+ 
  theme_light()


# PCA with correlation matrix
pca_gradesCor=prcomp(grades,scale.=T)
summary(pca_gradesCor)

# Loadings and variance of the components 
pca_gradesCor$rotation
pca_gradesCor$sdev^2

# Biplot
fviz_pca_biplot(pca_gradesCor,title="",
                ggtheme = theme_minimal(),geom="text")
###########################################################

###########################################################
# Example: Pizzas
pizza=read.csv(here("R/DimensionalityReduction/Data/Pizza.csv"),header = T)

# PCA with covariance matrix 
pca_pizza=prcomp(pizza[,-c(1,2)])
summary(pca_pizza)

# Loadings and variance of components
pca_pizza$rotation
pca_pizza$sdev^2

# Biplot
fviz_pca_biplot(pca_pizza,title="",
                habillage=pizza$brand,
                ggtheme = theme_minimal(),geom="point")

# PCA with correlation matrix
pca_pizzaCor=prcomp(pizza[,-c(1,2)],scale.=T)
summary(pca_pizzaCor)

# Loadings and variance of the components
pca_pizzaCor$rotation

x=c(1,2,3,4,5,6,7)
y=pca_pizzaCor$sdev^2;y
exp_var=data.frame("PC"=x,"Var"=y)

ggplot(data=exp_var,aes(x=PC,y=Var))+
  geom_col(fill="navyblue",alpha=.85)+ 
  theme_light()

# Biplot
fviz_pca_biplot(pca_pizzaCor,title="",
                habillage=pizza$brand,
                ggtheme = theme_minimal(),geom="point")
###########################################################