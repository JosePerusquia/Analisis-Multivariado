##################################################################
# Dirichlet distribution                               
# Author: Jose Antonio Perusquia Cortes
# Afil: Facultad de Ciencias - UNAM
# Module: Multivariate Analysis
##################################################################

##################################################################
# Required libraries                                                     
library(ggplot2)            # Version 4.0.2
library(GGally)             # Version 2.4.0
library(ggthemes)           # Version 5.2.0
library(scatterplot3d)      # Version 0.3-45
library(plotly)             # Version 4.12.0
library(here)               # Version 1.0.2
##################################################################

##################################################################
# Dirichlet Distribution Functions                                     

# Function that generates a random sample of size n of a Dirichlet
# Distribution with concentration parameter alpha
rDirichlet=function(n,alpha){
  p = length(alpha)
  Y=matrix(nrow=n,ncol=p)
  for(i in 1:p){
    Y[,i]=rgamma(n,alpha[i])
  }
  V = rowSums(Y)
  return(Y / V)
}

# Density of the Dirichlet distribution of concentration parameter
# alpha if log = T it returns the log-density.
dDirichlet=function(x,alpha,log=T){
  if(is.vector(x)) {
    x = matrix(x, nrow = 1)
  }
  n=dim(x)[1]
  res=numeric(0)
  
  c0 = lgamma(sum(alpha)) - sum(lgamma(alpha))
  
  logdens = apply(x,1,function(row)
    c0 + sum((alpha-1)*log(row))
  )
  
  if(log){
    return(logdens)
  }else{
    return(exp(logdens))
  }
}

# Function that plots a sample of a Dirichlet distribution of 
# dimension 3.
gDirichlet=function(x,fx){
  
  x_cart=x[,2]+.5*x[,3]
  y_cart=sqrt(3)*x[,3]/2
  
  df=data.frame(x=x_cart,y=y_cart,z=fx)
  
  p=ggplot(data=df,aes(x=x_cart,y=y_cart))+
    annotate('segment',x=0,y=0,xend=1,yend=0,linewidth=.5)+
    annotate('segment',x=0,y=0,xend=.5,yend=sqrt(3)/2,linewidth=.5)+
    annotate('segment',x=.5,y=sqrt(3)/2,xend=1,yend=0,linewidth=.5)+
    annotate("text", x = -.02, y = -.0025, label = list('X[1]'),
             parse = TRUE)+
    annotate("text", x = 1.02, y = -.0025, label = list('X[2]'),
             parse = TRUE)+
    annotate("text", x = .505, y = .89, label = list('X[3]'),
             parse = TRUE)+
    geom_point(aes(col=z),shape=1,size=.5,alpha=.7)+
    theme_void()+
    scale_colour_gradient(low='blue',high='red',
                          name=expression(f[x]))+
    labs(x='',y='')+
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          legend.position = c(.9,.75))
  plot(p)
}

# Heatmap of the whole distribution 
DirichletHeatmap=function(alpha,resolution=500){
  
  # Create grid
  x1=seq(0,1,length.out=resolution)
  x2=seq(0,1,length.out=resolution)
  
  grid=expand.grid(x1=x1,x2=x2)
  grid$x3=1-grid$x1-grid$x2
  
  # Keep only points inside the simplex
  grid=grid[grid$x3>=0,]
  
  X=as.matrix(grid[,c("x1","x2","x3")])
  
  # Evaluate density
  grid$dens=dDirichlet(X,alpha,log=FALSE)
  
  # Barycentric → Cartesian
  grid$x_cart=grid$x2 + 0.5*grid$x3
  grid$y_cart=sqrt(3)/2 * grid$x3
  
  # Plot
  p=ggplot(grid,aes(x=x_cart,y=y_cart,fill=dens))+
    geom_tile()+
    annotate('segment',x=0,y=0,xend=1,yend=0,linewidth=.5)+
    annotate('segment',x=0,y=0,xend=.5,yend=sqrt(3)/2,linewidth=.5)+
    annotate('segment',x=.5,y=sqrt(3)/2,xend=1,yend=0,linewidth=.5)+
    annotate("text",x=-.02,y=-.02,label="X[1]",parse=TRUE)+
    annotate("text",x=1.02,y=-.02,label="X[2]",parse=TRUE)+
    annotate("text",x=.5,y=.9,label="X[3]",parse=TRUE)+
    coord_equal()+
    theme_void()+
    scale_fill_viridis_c(name=expression(f[x]))+
    theme(legend.position=c(.9,.75))
  
  plot(p)
}

# A 3D plot using plotly to observe the density in the simplex
Dirichlet3D=function(alpha,resolution=120){
  
  # Create simplex grid
  x1=seq(0,1,length.out=resolution)
  x2=seq(0,1,length.out=resolution)
  
  grid=expand.grid(x1=x1,x2=x2)
  grid$x3=1-grid$x1-grid$x2
  
  # Keep simplex points
  grid=grid[grid$x3>=0,]
  
  X=as.matrix(grid[,c("x1","x2","x3")])
  
  # Evaluate density
  grid$dens=dDirichlet(X,alpha,log=FALSE)
  
  # Barycentric → Cartesian
  grid$x=grid$x2 + .5*grid$x3
  grid$y=sqrt(3)/2 * grid$x3
  grid$z=grid$dens
  
  # Plot interactive surface
  plot_ly(grid,
          x=~x,
          y=~y,
          z=~z,
          type="mesh3d",
          intensity=~z,
          colorscale="Viridis",
          showscale=TRUE) |>
    layout(scene=list(
      xaxis=list(title=""),
      yaxis=list(title=""),
      zaxis=list(title="Density"),
      aspectratio=list(x=1,y=1,z=0.7)
    ))
}
##################################################################

##################################################################
# Simulations and plots for different alphas
n=5000
k=3
alphas=matrix(c(1,1,1,5,5,5,1,2,2,2,4,8),byrow = T,nrow=4)

# All alphas set to 1 so it becomes uniform in the simplex
set.seed(314159)
x=rDirichlet(n,alphas[1,])
fx=dDirichlet(x,alphas[1,],log=F)
gDirichlet(x,fx)
DirichletHeatmap(alphas[1,])
Dirichlet3D(alphas[1,])

# All alphas set to 5 concentrates the sample in the middle
set.seed(314159)
x=rDirichlet(n,alphas[2,])
fx=dDirichlet(x,alphas[2,],log=F)
gDirichlet(x,fx)
DirichletHeatmap(alphas[2,])
Dirichlet3D(alphas[2,])

# alpha1 = 1 and alpha2=alpha3=2 concentrates the sample towards
# X2 and X3
set.seed(314159)
x=rDirichlet(n,alphas[3,])
fx=dDirichlet(x,alphas[3,],log=F)
gDirichlet(x,fx)
DirichletHeatmap(alphas[3,])
Dirichlet3D(alphas[3,])

# alpha1 = 2, alpha2 = 4 and alpha3 =8 concentrates the sample
# towards X3
set.seed(314159)
x=rDirichlet(n,alphas[4,])
fx=dDirichlet(x,alphas[4,],log=F)
gDirichlet(x,fx)
DirichletHeatmap(alphas[4,])
Dirichlet3D(alphas[4,])
#################################################################

