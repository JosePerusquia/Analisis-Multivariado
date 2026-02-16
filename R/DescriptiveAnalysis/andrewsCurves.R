############################################################################
# Título:
#     Curvas de Andrews
#
# Autor: 
#     José A. Perusquía Cortés
#
# Descripción: 
#     Código para generar las curvas de Andrews como fueron propuestas
#     en Andrews, D.F. (1972). Plots of High-Dimensional Data. Biometrics,
#     28(1), 125–136.
############################################################################

############################################################################
# Librerias
library(ggplot2)
############################################################################

andrewsCurves=function(x,groups,xlab="",ylab="",title="",
                       legend.title=element_blank(),
                       is.grouped=T){
  
  t=seq(-pi,pi,length.out=200)
  n=nrow(x)
  p=ncol(x)
  
  f=matrix(0,length(t),p)
  f[,1]=1/sqrt(2)
  
  for(i in 2:p){
    if(i%%2==0){
      f[,i]=sin(floor(i/2)*t)
    }else{
      f[,i]=cos(floor(i/2)*t)
    }
  }
  
  res=matrix(0,n,length(t))

  res[1,]=apply(t(f)*x[1,],2,sum)
  
  if(is.grouped){
    df=data.frame(x=t,y=res[1,],col=groups[1])
    
    p=ggplot(df)+
      geom_line(aes(x,y,colour=col))
    
    for(i in 2:n){
      res[i,]=apply(t(f)*x[i,],2,sum)
      df=data.frame(x=t,y=res[i,],col=groups[i])
      p=p+geom_line(data=df,mapping=aes(x,y,colour=col))
    }
    
    p=p+
      theme_minimal()+
      xlab(xlab)+
      ylab(ylab)+
      labs(title=title)+
      theme(
        legend.title = legend.title
      )
    plot(p)
    
  }else{
    df=data.frame(x=t,y=res[1,])
    
    p=ggplot(df)+
      geom_line(aes(x,y))
    
    for(i in 2:n){
      res[i,]=apply(t(f)*x[i,],2,sum)
      df=data.frame(x=t,y=res[i,])
      p=p+geom_line(data=df,mapping=aes(x,y))
    }
    
    p=p+
      theme_minimal()+
      xlab(xlab)+
      ylab(ylab)+
      labs(title=title)+
      theme(
        legend.position='none'
      )
    plot(p)
  }
}


