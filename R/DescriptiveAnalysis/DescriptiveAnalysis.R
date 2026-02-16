#################################################################
# Multivariate descriptive analysis                                           
# Author: Jose Antonio Perusquia Cortes
# Afil: Facultad de Ciencias - UNAM
# Module: Multivariate Analysis
#################################################################

#################################################################
# Required libraries                                                 
library(TeachingDemos)   # Version 2.13
library(pracma)          # Version 2.4.4
library(GGally)          # Version 2.2.1
library(corrplot)        # Version 0.95
library(here)            # Version 1.0.1
library(dplyr)           # Version 1.1.4
library(purrr)           # Version 1.0.4
#################################################################

#################################################################
# Source code for Andrews Curves using ggplot
source(here("R/DescriptiveAnalysis/andrewsCurves.R"))
#################################################################

#################################################################
# Iris data set

# Quick view on the characteristics of the data frame
str(iris)

# Mean
summary(iris)
apply(iris[,-5],2,mean)
colMeans(iris[,-5])
t(as.matrix(iris[,-5]))%*%rep(1,150)/150

# Mean by group
by(iris[,-5],iris[,5],colMeans)

# Covariance matrix
var(iris[,-5])
cov(iris[,-5])
W=as.matrix(sweep(iris[,-5],2,colMeans(iris[,-5])))
t(W)%*%W/(150-1)

# Covariance matrix by group
by(iris[,-5],iris[,5],cov)

# Correlation matrix
cor(iris[,-5])

# Correlation matrix by group
by(iris[,-5],iris[,5],cor)
#################################################################

#################################################################
# Dispersion plot using GGally
lower_fn = function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_point(alpha = 0.6) +
    theme_minimal()
}

diag_fn = function(data, mapping, ...) {
  ggally_densityDiag(data = data, mapping = mapping, ...) +
    theme_minimal()
}

ggpairs(
  iris[,-5],
  lower = list(continuous = lower_fn),
  diag = list(continuous = diag_fn)
)

# For grouped data using GGally
ggpairs(iris, aes(color = Species, alpha = 0.5),columns = 1:4,
        upper = list(continuous = wrap("cor", size = 2.5)),
        lower = list(continuous = lower_fn),
        diag = list(continuous = diag_fn))
#################################################################

#################################################################
# Correlation plot

# For a data set without labels
corrplot(cor(iris[,-5]),method="ellipse",
         tl.pos='l',tl.col='black')

# For grouped data we split the data and plot separately
iris %>%
  split(.$Species) %>%
  walk(~ corrplot(cor(select(.x, -Species)),
                  method = "ellipse",
                  tl.pos = "n"))
#################################################################

#################################################################
# Stars
stars(iris[,-5],lwd = 1)
#################################################################

#################################################################
# Chernoff's faces

# Original order of the variables
faces2(iris[,-5])

# The faces change if we change the order
faces2(iris[,c(4,3,2,1)])
#################################################################

#################################################################
# Andrews Curves

# For non-grouped data
andrewsCurves(as.matrix(iris[,-5]))

# For grouped data
andrewsCurves(as.matrix(iris[,-5]),iris[,5],
              legend.title = 'Species')

# The order changes the curves
andrewsCurves(as.matrix(iris[,c(4,3,2,1)]),iris[,5],
              legend.title = 'Species')

# Usin pracma library 
andrewsplot(as.matrix(iris[,-5]),iris[,5])
andrewsplot(as.matrix(iris[,-5]),iris[,5],style="cart")
#################################################################