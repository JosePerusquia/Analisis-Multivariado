################################################################################
# Libraries                                                                    #
################################################################################
library(ggplot2)
library(GGally)
library(ggthemes)
library(psych)
library(here)
library(expm)
library(GPArotation)

################################################################################
# Calificaciones Ejemplo                                                       

# Leer los datos
calificaciones=read.table(here("../Datos/calificaciones.txt"),header = T)

# Obtenemos matriz de correlación muestral
R=round(cor(calificaciones),3);R

# Realizamos factor de análisis principales sin rotación con 1 factor
efa=fa(r=R,nfactors=1,fm="pa",rotate="none")
efa

# Realizamos factor de análisis principales sin rotación con 2 facores
efa=fa(r=R,nfactors=2,fm="pa",rotate="none")
efa

# Realizamos factor de análisis por mle sin rotación con 1 factor
efa=fa(r=R,nfactors=1,n.obs=88,fm="mle",rotate="none")
efa

# Estadístico de prueba para 1 factor
efa$STATISTIC

# De forma alternativa construimos el estadístico de prueba para un factor
n=88;p=5;k=1

Lambda=efa$loadings
Psi=diag(efa$uniquenesses)
Sigma=(Lambda%*%t(Lambda))+Psi

minF=tr(solve(Sigma)%*%R)-log(det(solve(Sigma)%*%R))-p

n_prime=n-1-((2*p+5)/6)-((2*k)/3)

U=n_prime*minF;U

#Valor crítico
nu=.5*(p-k)^2-.5*(p+k)
qchisq(.95,nu)

# Realizamos factor de análisis por mle sin rotación con 2 factores
efa=fa(r=R,nfactors=2,n.obs=88,fm="mle",rotate="none")
efa

# Estadístico de prueba para 2 factores
efa$STATISTIC

# De forma alternativa construimos el estadístico de prueba para un factor
n=88;p=5;k=2

Lambda=efa$loadings
Psi=diag(efa$uniquenesses)
Sigma=(Lambda%*%t(Lambda))+Psi

minF=tr(solve(Sigma)%*%R)-log(det(solve(Sigma)%*%R))-p

n_prime=n-1-((2*p+5)/6)-((2*k)/3)

U=n_prime*minF;U

#Valor crítico
nu=.5*(p-k)^2-.5*(p+k)
qchisq(.95,1)

# Realizamos factor de análisis principales con rotación varimax
efa1=fa(r=R,nfactors=2,fm="pa",rotate="varimax")
efa1

# Graficamos la rotación
efa_df=data.frame("F1"=efa$loadings[,1],"F2"=efa$loadings[,2],
                  "Col"=as.factor(c("A","A","C","C","C")))

p=ggplot(efa_df,aes(x=F1,y=F2,colour=Col))+
  geom_point(size=1,alpha=0)+
  geom_text(show.legend=F,data=efa_df,
            aes(label = c("Lin.","Est.","Prob.","Fin.","Cal.")))+
  theme_minimal()+
  guides(colour = guide_legend(title="Examen", override.aes=list(alpha=1)))
  
p

# Primer eje

x_seq=seq(0,1,by=.1)
y1_seq=numeric(length(x_seq))
y2_seq=numeric(length(x_seq))

for(i in 1:length(x_seq)){
  y1_seq[i]=x_seq[i]*efa1$rot.mat[1,1]
  y2_seq[i]=x_seq[i]*efa1$rot.mat[2,1]
}

df_seq_x=data.frame("x1"=y1_seq,"x2"=y2_seq)
p=p+geom_line(data=df_seq_x,aes(x=x1,y=x2),col="red")
p

# Segundo eje

for(i in 1:length(x_seq)){
  y1_seq[i]=x_seq[i]*efa1$rot.mat[1,2]
  y2_seq[i]=x_seq[i]*efa1$rot.mat[2,2]
}

df_seq_x=data.frame("x1"=y1_seq,"x2"=y2_seq)
p=p+geom_line(data=df_seq_x,aes(x=x1,y=x2),col="red")
p

# Factores con rotación
efa_df=data.frame("F1"=efa1$loadings[,1],"F2"=efa1$loadings[,2],
                  "Col"=as.factor(c("A","A","C","C","C")))

p=ggplot(efa_df,aes(x=F1,y=F2,colour=Col))+
  geom_point(size=1,alpha=0)+
  geom_text(show.legend=F,data=efa_df,
            aes(label = c("Lin.","Est.","Prob.","Fin.","Cal.")))+
  theme_minimal()+
  guides(colour = guide_legend(title="Examen", override.aes=list(alpha=1)))

p


# Realizamos factor de análisis principales con rotación oblimin
efa1=fa(r=R,nfactors=2,fm="pa",rotate="oblimin")
efa1

# Graficamos la rotación
efa_df=data.frame("F1"=efa$loadings[,1],"F2"=efa$loadings[,2],
                  "Col"=as.factor(c("A","A","C","C","C")))

p=ggplot(efa_df,aes(x=F1,y=F2,colour=Col))+
  geom_point(size=1,alpha=0)+
  geom_text(show.legend=F,data=efa_df,
            aes(label = c("Lin.","Est.","Prob.","Fin.","Cal.")))+
  theme_minimal()+
  guides(colour = guide_legend(title="Examen", override.aes=list(alpha=1)))

p

# Primer eje
x_seq=seq(0,1,by=.1)
y1_seq=numeric(length(x_seq))
y2_seq=numeric(length(x_seq))

for(i in 1:length(x_seq)){
  y1_seq[i]=x_seq[i]*efa1$rot.mat[1,1]
  y2_seq[i]=x_seq[i]*efa1$rot.mat[2,1]
}

df_seq_x=data.frame("x1"=y1_seq,"x2"=y2_seq)
p=p+geom_line(data=df_seq_x,aes(x=x1,y=x2),col="red")
p

#Segundo eje
for(i in 1:length(x_seq)){
  y1_seq[i]=x_seq[i]*efa1$rot.mat[1,2]
  y2_seq[i]=x_seq[i]*efa1$rot.mat[2,2]
}

df_seq_x=data.frame("x1"=y1_seq,"x2"=y2_seq)
p=p+geom_line(data=df_seq_x,aes(x=x1,y=x2),col="red")
p

#Factores con rotación
efa_df=data.frame("F1"=efa1$loadings[,1],"F2"=efa1$loadings[,2],
                  "Col"=as.factor(c("A","A","C","C","C")))

p=ggplot(efa_df,aes(x=F1,y=F2,colour=Col))+
  geom_point(size=1,alpha=0)+
  geom_text(show.legend=F,data=efa_df,
            aes(label = c("Lin.","Est.","Prob.","Fin.","Cal.")))+
  theme_minimal()+
  guides(colour = guide_legend(title="Examen", override.aes=list(alpha=1)))

p
################################################################################
