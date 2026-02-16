###################################################################
# Source code for the Andrews curves as proposed by Andrews (1972)
# Author: Jose Antonio Perusquia Cortes
# Afil : Facultad de Ciencias - UNAM
# Module : Multivariate Analysis
###################################################################

###################################################################
# Required libraries
library(ggplot2)          # Version 3.5.2
library(ggthemes)         # Version 5.1.0
###################################################################

###################################################################
# Source code for the original implementation of the Andrews 
# Curves. The function receives as parameters:
# x = the matrix of observations
# groups = the labels for each observation
# xlab, ylab, title, legend.title for the plot

andrewsCurves <- function(x, groups = NULL,
                          xlab = "", ylab = "", title = "",
                          legend.title = NULL) {
  
  x <- as.matrix(x)
  
  t <- seq(-pi, pi, length.out = 200)
  n <- nrow(x)
  p <- ncol(x)
  
  # Build Fourier basis
  f <- matrix(0, length(t), p)
  f[,1] <- 1/sqrt(2)
  
  for(i in 2:p){
    if(i %% 2 == 0){
      f[,i] <- sin((i/2)*t)
    } else {
      f[,i] <- cos(((i-1)/2)*t)
    }
  }
  
  # Vectorized transformation
  res <- x %*% t(f)
  
  # Long format
  df <- data.frame(
    t = rep(t, each = n),
    value = as.vector(res),
    id = rep(1:n, times = length(t))
  )
  
  if(!is.null(groups)){
    df$group <- rep(groups, times = length(t))
  }
  
  # Base plot
  if(!is.null(groups)){
    p <- ggplot(df, aes(t, value, group = id, colour = group)) +
      geom_line()
    
    # Apply legend title if provided
    if(!is.null(legend.title)){
      p <- p + labs(colour = legend.title)
    }
    
  } else {
    p <- ggplot(df, aes(t, value, group = id)) +
      geom_line() +
      theme(legend.position = "none")
  }
  
  p <- p +
    theme_minimal() +
    labs(x = xlab, y = ylab, title = title)
  
  return(p)
}
###################################################################