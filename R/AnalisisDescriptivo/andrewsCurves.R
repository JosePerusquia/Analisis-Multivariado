library(ggplot2)

andrewsCurves=function(x,groups){
  
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
    xlab("")+
    ylab("")+
    theme(
      legend.title = element_blank()
    )
  plot(p)
  
}


