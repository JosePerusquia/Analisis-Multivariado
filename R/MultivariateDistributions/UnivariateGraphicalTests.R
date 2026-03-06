############################################################################
# Univariate graphical normality tests
# Author: José A. Perusquía Cortés
# Afil: Facultad de Ciencias - UNAM
# Module: Multivariate Analysis
############################################################################

############################################################################
# Required libraries
library(ggplot2)          # Version 4.0.2
library(ggpubr)           # Version 0.6.3
library(ggthemes)         # Version 5.2.0
library(nortest)          # Version 1.0-4
############################################################################

##############################################################################
# Function that plots histogram, qqplot, boxplot and one of the 
# nonparametric normality tests supported by nortest library. By 
# default it performs Anderson-Darling.
univariateNormalityPlots = function(X, test = "ad") {
  
  library(ggplot2)
  library(patchwork)
  library(nortest)
  
  n <- nrow(X)
  p <- ncol(X)
  var_names <- colnames(X)
  
  for(i in 1:p){
    
    x <- X[,i]
    df <- data.frame(x = x)
    
    ## Histogram with normal curve
    h <- ggplot(df, aes(x)) +
      geom_histogram(aes(y = after_stat(density)),
                     bins = 15,
                     fill = "skyblue3",
                     color = "black",
                     alpha = .6) +
      stat_function(fun = dnorm,
                    args = list(mean = mean(x),
                                sd = sd(x)),
                    color = "darkred",
                    linewidth = 1) +
      theme_minimal() +
      labs(title = paste("Histogram:", var_names[i]),
           x = "", y = "")
    
    
    ## QQ plot
    q <- ggplot(df, aes(sample = x)) +
      stat_qq(size = 1) +
      stat_qq_line(color = "red") +
      theme_minimal() +
      labs(title = "Q-Q Plot", x = "", y = "")
    
    
    ## Boxplot
    bx <- ggplot(df, aes(y = x)) +
      geom_boxplot(fill = "skyblue3",
                   color = "black",
                   outlier.alpha = .7) +
      theme_minimal() +
      labs(title = "Boxplot", x = "", y = "")
    
    
    ## Normality test
    p_val <- switch(test,
                    ad = nortest::ad.test(x)$p.value,
                    cvm = nortest::cvm.test(x)$p.value,
                    lillie = nortest::lillie.test(x)$p.value,
                    sf = nortest::sf.test(x)$p.value)
    
    label <- paste0(toupper(test),
                    " test p-value = ",
                    round(p_val,5))
    
    pv <- ggplot() +
      annotate("text",
               x = 0,
               y = 0,
               label = label,
               size = 5) +
      theme_void()
    
    
    ## Combine plots
    figure <- (h | q) /
      (bx | pv)
    
    print(figure)
    
  }
}
################################################################################

################################################################################
# Function that plots the qqplot of the quadratic form of a normal 
# distribution against the quantiles of a chi-squared distribution of 
# p degrees of freedom. By default, the function estimates the mean
# and the covariance matrix, otherwise they need to be supplied
quadraticFormPlot = function(X, estimate = TRUE, mu = NULL, Sigma = NULL){
  
  library(ggplot2)
  
  n <- nrow(X)
  p <- ncol(X)
  
  if(estimate){
    mu <- colMeans(X)
    Sigma <- cov(X)
  }
  
  ## Mahalanobis distances
  D2 <- mahalanobis(X, mu, Sigma)
  
  ## theoretical quantiles
  probs <- ppoints(n)
  theo <- qchisq(probs, df = p)
  
  ## empirical quantiles
  emp <- sort(D2)
  
  df <- data.frame(theo, emp)
  
  q <- ggplot(df, aes(x = theo, y = emp)) +
    geom_point(size = 1) +
    geom_abline(slope = 1, intercept = 0, color = "red") +
    theme_minimal() +
    labs(x = expression(chi[p]^2~"theoretical quantiles"),
         y = "Empirical Mahalanobis distances",
         title = "Chi-square Q-Q Plot")
  
  print(q)
}
################################################################################


