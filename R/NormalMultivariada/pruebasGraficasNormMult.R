############################################################################
# Título: Pruebas normalidad multivariada
# Autor: José A. Perusquía Cortés
# Afil: Facultad de Ciencias - UNAM
# Curso: Análisis Multivariado
#
# Descripción: Código para generar las gráficas básicas para checar
# la normalidad multivariada de un conjunto de datos
############################################################################

############################################################################
# Librerias
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(nortest)
############################################################################

##############################################################################
# Función para crear el histograma, el qqplot, boxplot y
# alguna de las pruebas no paramétricas soportadas por la librería nortest
# (exceptuando la de Pearson). Por default realiza la de Anderson-Darling 'ad'
# y se puede cambiar por Cramer-von Mises 'cvm', Lilliefors 'lillie' o
# Shapiro-Francia 'sf'.
pruebas_normalidad_univariadas=function(X,prueba="ad"){
  n=nrow(X)
  p=ncol(X)
  names=colnames(X)
  
  for(i in 1:p){
    df=data.frame(var=X[,i])
    
    # Histograma
    h=ggplot(data=df,aes(var))+
      geom_histogram(breaks=hist(X[,i],plot=F)$breaks,alpha=.5,
                     col='black',fill='skyblue3')+
      labs(title='',x='',y='')+
      theme_minimal()+
      theme(axis.title =element_text(size=8))
    
    # qqnorm  
    probs=ppoints(n,.5)
    emp_quantiles=quantile(X[,i],probs=probs)
    the_quantiles=qnorm(probs)
    
    x1=qnorm(.25)
    x2=qnorm(.75)
    y1=quantile(X[,i],.25)
    y2=quantile(X[,i],.75)
    
    m=(y2-y1)/(x2-x1)
    b=y1-m*x1
    
    quantiles=data.frame(emp_quantiles,the_quantiles)
    names(quantiles)=c('Emp','Teo')
    
    q=ggplot(quantiles,aes(x=Teo,y=Emp))+
      geom_point(col='black',size=1)+
      geom_abline(intercept = b ,slope =m,col="red")+
      labs(title='',x='',y='')+
      theme_minimal()+
      theme(axis.title =element_text(size=8))
    
    # boxplot
    bx=ggplot(df, aes(var)) +
      geom_boxplot(outlier.alpha=.75,col='black',fill='skyblue3')+
      theme_minimal()+
      labs(x="")+
      ylim(-.75,.75)
    
    #pruebas normalidad
    if(prueba=="ad"){
      p_val=ad.test(X[,i])$p
      prueba_res=paste("Anderson-Darling:",' ',round(p_val,5))
      x=0
      y=0
      df1=data.frame(x=x,y=y,text=prueba_res)
      pv=ggplot(df1, aes(x,y,label=text)) +
        geom_text()+
        theme_transparent()
    }else if(prueba=="cvm"){
      p_val=cvm.test(X[,i])$p
      prueba_res=paste("Cramer-von Mises:",' ',round(p_val,5))
      x=0
      y=0
      df1=data.frame(x=x,y=y,text=prueba_res)
      pv=ggplot(df1, aes(x,y,label=text)) +
        geom_text()+
        theme_transparent()
    }else if(prueba=="lillie"){
      p_val=lillie.test(X[,i])$p
      prueba_res=paste("Lilliefors:",' ',round(p_val,5))
      x=0
      y=0
      df1=data.frame(x=x,y=y,text=prueba_res)
      pv=ggplot(df1, aes(x,y,label=text)) +
        geom_text()+
        theme_transparent()
    }else if(prueba=="sf"){
      p_val=sf.test(X[,i])$p
      prueba_res=paste("Shapiro-Francia:",' ',round(p_val,5))
      x=0
      y=0
      df1=data.frame(x=x,y=y,text=prueba_res)
      pv=ggplot(df1, aes(x,y,label=text)) +
        geom_text()+
        theme_transparent()
    }

    figure = ggarrange(h,q,bx,pv,
                       labels = c("", "","",""),
                       ncol = 2, nrow = 2)

    plot(figure)
      
  }
  

  
}
################################################################################

################################################################################
# Función para crear el qqplot de la forma cuadrática de la distribución
# normal multivariada para compararla con los cuantiles de una distribución
# ji-cuadrada de p grados de libertad. Por default la función hace estimación
# del vector de medias y de la matriz de covarianzas, de lo contrario se
# requiere proveer el vector de medias y la matriz de covarianzas teórica
prueba_forma_cuad=function(X,estimate="T",mu,sigma){
  
  n=nrow(X)
  p=ncol(X)
  
  chi_sq=numeric(n)
  
  if(estimate=="T"){
    mu=colMeans(X)
    sigma=solve(cov(X))
  }
  
  for(i in 1:n){
    chi_sq[i]=t(X[i,]-mu)%*%sigma%*%(X[i,]-mu)
  }
  
  probs=ppoints(n,.5)
  
  emp_quantiles=quantile(chi_sq,probs=ppoints(n,.5))
  the_quantiles=qchisq(ppoints(n,.5),p)
  
  x1=qchisq(.25,p)
  x2=qchisq(.75,p)
  
  y1=quantile(chi_sq,.25)
  y2=quantile(chi_sq,.75)
  
  m=(y2-y1)/(x2-x1)
  b=y1-m*x1
  
  quantiles=data.frame(emp_quantiles,the_quantiles)
  names(quantiles)=c('Teo','Emp')
  
  q=ggplot(quantiles,aes(x=Emp,y=Teo))+
    geom_point(col='black',size=1)+
    geom_abline(intercept = b ,slope =m,col="red")+
    labs(title='',x='',y='')+
    theme_minimal()+
    theme(axis.title =element_text(size=8))
  plot(q)
}
################################################################################


