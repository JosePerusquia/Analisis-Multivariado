##################################################################
# Multivariate t and multivariate Cauchy distribution                               
# Author: Jose Antonio Perusquia Cortes
# Afil: Facultad de Ciencias - UNAM
# Module: Multivariate Analysis
##################################################################

##################################################################
# Required libraries                                                     
library(ggplot2)            # Version 4.0.2
library(GGally)             # Version 2.4.0
library(ggthemes)           # Version 5.2.0
library(mvtnorm)            # Version 1.3-3
library(scatterplot3d)      # Version 0.3-45
library(expm)               # Version 1.0-0
library(plotly)             # Version 4.12.0
library(MVN)                # Version 6.3
library(here)               # Version 1.0.2
##################################################################

##################################################################
# Source code for functions used to perform a preliminar analysis
# on the distributional assumption of multivariate normality
source(here('UnivariateGraphicalTests.R'))
##################################################################

##################################################################
# Multivariate t distribution centred in (0,0) and with covariance
# matrix given by sigma
delta=c(0,0)
sigma=matrix(c(3,1,1,3),byrow = T,nrow=2)

# Density
x = seq(-6, 6 , length.out = 100)
y = seq(-6 ,6, length.out = 100)
grid = expand.grid(x=x,y=y)
mvds = dmvt(grid,sigma=sigma,df=2,log=F)
z_matrix = matrix(mvds, nrow = length(x), byrow = FALSE)

persp(x, y, z_matrix,
      theta = 180, phi = 30,
      expand = .5, shade = 0.1,
      col = "lightblue",
      xlab = expression(x),
      ylab = expression(y),
      zlab = expression(f))

plot_ly(x=x,y=y,z=z_matrix,type="surface")

# Random sample and density for scatterplot
multt.sample=rmvt(1000,sigma,df=2)
dens=dmvt(multt.sample,sigma=sigma,df=2,log=F)
scatterplot3d(multt.sample[,1],
              multt.sample[,2],
              dens,pch=1,
              color="blue",cex.symbols = 0.7,
              type="p",xlab=expression(x[1]),
              ylab=expression(x[2]),zlab=expression(f(x)))

# Contour levels with ggplot
grid$z = mvds
ggplot(grid, aes(x, y, z = z)) +
  geom_contour(col='cyan4') +
  coord_fixed() +
  theme_minimal() +
  labs(x='',y='')

# Contour levels with plotly
plot_ly(x = x,y = y,z = z_matrix,type = "contour")

# Mardia test
set.seed(31415)
multt.sample=rmvt(50,sigma,df=5)
univariateNormalityPlots(multt.sample)
quadraticFormPlot(multt.sample)

mvn_HZ=mvn(multt.sample,mvn_test="hz")
mvn_HZ$multivariate_normality

mvn_RS=mvn(multt.sample,mvn_test="royston")
mvn_RS$multivariate_normality

mvn_MA=mvn(multt.sample,mvn_test="mardia")
mvn_MA$multivariate_normality
##################################################################

##################################################################
# Multivariate Cauchy distribution centred in (0,0) and with 
# covariance matrix given by sigma
delta=c(0,0)
sigma=matrix(c(3,1,1,3),byrow = T,nrow=2)

# Density
x = seq(-6, 6 , length.out = 100)
y = seq(-6 ,6, length.out = 100)
grid = expand.grid(x=x,y=y)
mvds = dmvt(grid,sigma=sigma,df=1,log=F)
z_matrix = matrix(mvds, nrow = length(x), byrow = FALSE)

persp(x, y, z_matrix,
      theta = 180, phi = 30,
      expand = .5, shade = 0.1,
      col = "lightblue",
      xlab = expression(x),
      ylab = expression(y),
      zlab = expression(f))

plot_ly(x=x,y=y,z=z_matrix,type="surface")

# Random sample and density for scatterplot
multt.sample=rmvt(1000,sigma,df=1)
dens=dmvt(multt.sample,sigma=sigma,df=1,log=F)
scatterplot3d(multt.sample[,1],
              multt.sample[,2],
              dens,pch=1,
              color="blue",cex.symbols = 0.7,
              type="p",xlab=expression(x[1]),
              ylab=expression(x[2]),zlab=expression(f(x)))

# Contour levels with ggplot
grid$z = mvds
ggplot(grid, aes(x, y, z = z)) +
  geom_contour(col='cyan4') +
  coord_fixed() +
  theme_minimal() +
  labs(x='',y='')

# Contour levels with plotly
plot_ly(x = x,y = y,z = z_matrix,type = "contour")

# Mardia test
set.seed(31415)
multt.sample=rmvt(50,sigma,df=1)
univariateNormalityPlots(multt.sample)
quadraticFormPlot(multt.sample)

mvn_HZ=mvn(multt.sample,mvn_test="hz")
mvn_HZ$multivariate_normality

mvn_RS=mvn(multt.sample,mvn_test="royston")
mvn_RS$multivariate_normality

mvn_MA=mvn(multt.sample,mvn_test="mardia")
mvn_MA$multivariate_normality
##################################################################