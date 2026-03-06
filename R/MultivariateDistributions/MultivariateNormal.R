##################################################################
# Multivariate normal distribution                                
# Author: Jose Antonio Perusquia Cortes
# Afil: Facultad de Ciencias - UNAM
# Module: Multivariate Analysis
##################################################################

##################################################################
# Required libraries   
library(ggplot2)            # Version 3.5.2
library(GGally)             # Version 2.2.1
library(ggthemes)           # Version 5.1.0
library(mvtnorm)            # Version 1.3-3
library(scatterplot3d)      # Version 0.3-44
library(expm)               # Version 1.0-0
library(plotly)             # Version 4.11.0
library(MVN)                # Version 6.3
library(car)                # Version 3.1-3
library(ICSNP)              # Version 1.1-2
library(here)               # Version 1.0.1
##################################################################

##################################################################
# Source code for functions used to perform a preliminar analysis
# on the distributional assumption of multivariate normality
source(here('UnivariateGraphicalTests.R'))
##################################################################

##################################################################
# Non-degenerate bivariate normal distribution                                               
mu=c(0,0)
sigma=matrix(c(3,1,1,3),byrow = T,nrow=2)

# Density
x = seq(-6, 6 , length.out = 100)
y = seq(-6 ,6, length.out = 100)
grid = expand.grid(x=x,y=y)
mvds = dmvnorm(grid,mu,sigma)
z_matrix = matrix(mvds, nrow = length(x), byrow = FALSE)

persp(x, y, z_matrix,
      theta = 180, phi = 30,
      expand = .5, shade = 0.1,
      col = "lightblue",
      xlab = expression(x),
      ylab = expression(y),
      zlab = expression(f))

# Random sample and density for scatterplot
multnorm.sample=rmvnorm(1000,mu,sigma)
dens=dmvnorm(multnorm.sample,mean=mu,sigma=sigma)
scatterplot3d(multnorm.sample[,1],
              multnorm.sample[,2],
              dens,pch=1,
              color="blue",cex.symbols = 0.7,
              type="p",xlab=expression(x[1]),
              ylab=expression(x[2]),zlab=expression(f(x)))

# Transformation to a standard multivariate normal distribution
X=rmvnorm(10000,mu,sigma)
colnames(X) = c("X1","X2")
ggpairs(as.data.frame(X),
        upper = list(continuous = wrap("cor",size = 9)))

A = sqrtm(solve(sigma))
Y = t(A %*% t(X - matrix(mu, nrow(X), 2, byrow = TRUE)))
colnames(Y) = c('Y1','Y2')
ggpairs(as.data.frame(Y),
        upper = list(continuous = wrap("cor",size=9)))

# Contour levels with ggplot
grid$z = mvds
ggplot(grid, aes(x, y, z = z)) +
  geom_contour(col='cyan4') +
  coord_fixed() +
  theme_minimal() +
  labs(x='',y='')

# Contour levels with plotly
plot_ly(x = x,y = y,z = z_matrix,type = "contour")
##################################################################

##################################################################
# Examples on multivariate normality

# Simulated data univariate and multivariate tests
univariateNormalityPlots(multnorm.sample)
quadraticFormPlot(multnorm.sample)

mvn_HZ=mvn(multnorm.sample,mvn_test="hz")
mvn_HZ$multivariate_normality

mvn_RS=mvn(multnorm.sample,mvn_test="royston")
mvn_RS$multivariate_normality

mvn_MA=mvn(multnorm.sample,mvn_test="mardia")
mvn_MA$multivariate_normality

# Iris univariate and multivariate normality tests
univariateNormalityPlots(as.matrix(iris[,-5]))
quadraticFormPlot(as.matrix(iris[,-5]))

mvn_HZ = mvn(as.matrix(iris[,-5]),mvn_test="hz")
mvn_HZ$multivariate_normality

mvn_RS = mvn(as.matrix(iris[,-5]),mvn_test="royston")
mvn_RS$multivariate_normality

mvn_MA = mvn(as.matrix(iris[,-5]),mvn_test="mardia")
mvn_MA$multivariate_normality
##################################################################

##################################################################
# Randnomness of Wishart distribution through the ellipses             
set.seed(3141592)
W=rWishart(n=4,df=2,Sigma=diag(1,2))

ellipse(c(0, 0), shape=W[,,1], radius=1, col="red", 
        lty=2,add=F,xlim=c(-2,2),ylim=c(-2,2))
ellipse(c(0, 0), shape=W[,,2], radius=1, col="blue", lty=2)
ellipse(c(0, 0), shape=W[,,3], radius=1, col="green", lty=2)
ellipse(c(0, 0), shape=W[,,4], radius=1, col="purple", lty=2)
##################################################################

##################################################################
# Hypothesis testing for the mean of one population                                   

# Simulate the data
set.seed(32131)
mu=c(64.1,64.7)
sigma=matrix(c(191,155.6,155.6,313.5),byrow = T,nrow=2)
multnorm.sample=data.frame(rmvnorm(203,mu,sigma))
n = nrow(multnorm.sample)
p = ncol(multnorm.sample)

# Test Ho:mu = (60,60) vs Ha: mu != (60,60)
mu_0 = c(60,60)

# Known Sigma
p1=ggplot(multnorm.sample,aes(x=X1,y=X2))+
  geom_point()+
  geom_point(x=mu_0[1],y=mu_0[2],col="red",shape=3)+
  labs(x=expression(X[1]),y=expression(X[2]))+
  theme_light()
p1

# Sample mean
mu_bar=colMeans(multnorm.sample);mu_bar

# Statistic
T2_known = n*t(mu_bar-mu_0)%*%solve(sigma)%*%(mu_bar-mu_0)
T2_known

# Critical value at alpha=.05
qchisq(.95,p)

# Confidence region
ConfReg=data.frame(ellipse(mu_bar, shape=sigma, 
                           radius=sqrt(qchisq(.95,p)/n), 
                           col="red",add=F,draw=F, lty=2))

p1+
  geom_point(data=ConfReg,aes(x=x,y=y),col="blue",size=.1)+
  geom_point(x=mu_bar[1],y=mu_bar[2],col='blue',shape=3)
  
# Unknown Sigma
sigma_hat=cov(multnorm.sample);sigma_hat

# Statistic
T2_unknown=((n-p)*n)/(p*(n-1))*t(mu_bar-mu_0)%*%solve(sigma_hat)%*%(mu_bar-mu_0)
T2_unknown

# Critical value at alpha=.05
qf(.95,p,n-p)

# Confidence region
ConfReg=data.frame(ellipse(mu_bar,shape=sigma_hat, 
                    radius=sqrt((p*(n-1))/((n-p)*n)*qf(.95,p,n-p)), 
                    col="red",add=F,draw=F, lty=2))

p1+
  geom_point(data=ConfReg,aes(x=x,y=y),col="blue",size=.1)+
  geom_point(x=mu_bar[1],y=mu_bar[2],col='blue',shape=3)

# In R the ICSNP library contains the function HotellingsT2 that 
# performs one and two-sample procedures
HotellingsT2(multnorm.sample,mu=mu_0)
##################################################################

##################################################################
# Hypothesis testing for comparing the mean of two independent
# normal populations

# Iris virginica species and multivariate normality tests
irisVirginica = iris%>%
  filter(Species=='virginica')%>%
  select('Sepal.Length','Sepal.Width')

mvn_HZ = mvn(as.matrix(irisVirginica),mvn_test="hz")
mvn_HZ$multivariate_normality

mvn_RS = mvn(as.matrix(irisVirginica),mvn_test="royston")
mvn_RS$multivariate_normality

mvn_MA = mvn(as.matrix(irisVirginica),mvn_test="mardia")
mvn_MA$multivariate_normality

# Iris versicolor species and multivariate normality tests
irisVersicolor = iris%>%
  filter(Species=='versicolor')%>%
  select('Sepal.Length','Sepal.Width')

mvn_HZ = mvn(as.matrix(irisVersicolor),mvn_test="hz")
mvn_HZ$multivariate_normality

mvn_RS = mvn(as.matrix(irisVersicolor),mvn_test="royston")
mvn_RS$multivariate_normality

mvn_MA = mvn(as.matrix(irisVersicolor),mvn_test="mardia")
mvn_MA$multivariate_normality

# Hypothesis testing on the equality of covariance matrices
n = nrow(irisVersicolor)
m = nrow(irisVirginica)
p = 2

S1 = var(irisVersicolor);S1
Q1 = (n-1)*S1

S2 = var(irisVirginica);S2
Q2 = (m-1)*S2

# Statistic
t1 = (n+m)^((n+m)*p/2)/((n^(n*p/2))*(m^(m*p/2)))
t2 = det(Q1)^(n/2)*det(Q2)^(m/2)/(det(Q1+Q2)^((n+m)/2))
Q = -2*log(t1*t2);Q

# Compare it againts chi-squared distribution with p(p*1)/2 df
v = p*(p+1)/2
qchisq(.95,v)

# Since we do not reject the null hypothesis on the equality of 
# covariance matrices we perform a two sample T2 Hotelling test
Su = ((n-1)*S1 + (m-1)*S2)/(n+m-2)

xbar = colMeans(irisVersicolor)
ybar = colMeans(irisVirginica)

# Statistic
c1 = (n+m-p-1)*(n*m)/(p*(n+m-2)*(n+m))
c2 = t(ybar-xbar)%*%solve(Su)%*%(ybar-xbar)

delta2 = c1*c2;delta2

# Compare it against F distribution with p,n+m-p-1 df
qf(.95,p,n+m-p-1)

# Alternatively HotellingsT2 function
HotellingsT2(irisVersicolor,irisVirginica)
##################################################################
