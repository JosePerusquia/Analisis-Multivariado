################################################################################
# Análisis de correspondencias                                        
# Autor: Jose Antonio Perusquia Cortes
# Afil : Facultad de Ciencias - UNAM
# Curso : Análisis Multivariado
################################################################################

##############################################################################
# Librerías
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(expm)
library(ca)
library(here)
library(dplyr)
library(readxl)
##############################################################################

##############################################################################
# Datos estado de salud
health=matrix(c(243,789,167,18,6,
                220,809,164,35,6,
                147,658,181,41,8,
                90,469,236,50,16,
                53,414,306,106,30,
                44,267,284,98,20,
                20,136,157,66,17),
              nrow=7,ncol=5,byrow=T)

## Total de observaciones
n=sum(health);n

## Suma por renglones y columnas
row_sum=apply(health,1,sum);row_sum
col_sum=apply(health,2,sum);col_sum

## Masas y centroide
row_masses=row_sum/n;row_masses
centroide=col_sum/n;centroide

## Matrices diagonales de masas y centroide
Dr=diag(row_masses);Dr
Dc=diag(1/centroide);Dc

## Raiz cuadrada de matrices diagonales
Dr_sq=sqrtm(Dr);Dr_sq
Dc_sq=sqrtm(Dc);Dc_sq

## Perfiles por renglón
profiles=matrix(0,nrow=7,ncol=5)
for(i in 1:7){
  profiles[i,]=health[i,]/row_sum[i]
}

profiles

## Matriz quitando el centroide
R=matrix(0,nrow=7,ncol=5)
for(i in 1:7){
  R[i,]=profiles[i,]-centroide
}
R

## Encontramos la descomposición gsvd 
res=svd(Dr_sq%*%R%*%Dc_sq)

## Valores propios
values=res$d;values

## Matrices
U=res$u;U
V=res$v;V

N=solve(Dr_sq)%*%U;N
M=solve(Dc_sq)%*%V;M

## Coordenadas de los renglones
f=N[,c(1:2)]%*%diag(values[c(1:2)]);f

df_points=data.frame(x=f[,1],y=f[,2])
df_labels=data.frame(x=f[,1],y=f[,2])

ggplot(data=df_points,aes(x=x,y=y))+
  geom_text(show.legend=F,size=3,
            aes(label = c("16-24","25-34","35-44",
                          "45-54","55-64","65-74","75+")))+
  theme_minimal()+
  labs(x="",y="")+
  ylim(c(-.2,.2))
##############################################################################

##############################################################################
# Problema dual estado de salud
health_d=t(health)

## Suma por renglones y columnas
row_sum_d=apply(health_d,1,sum)
col_sum_d=apply(health_d,2,sum)

## Masas y centroide
row_masses_d=row_sum_d/n
centroide_d=col_sum_d/n

## Matrices diagonales de masas y centroide
Dr_d=diag(row_masses_d)
Dc_d=diag(1/centroide_d)

## Raiz cuadrada de matrices diagonales
Dr_sq_d=sqrtm(Dr_d)
Dc_sq_d=sqrtm(Dc_d)

## Perfiles por renglón
profiles_d=matrix(0,nrow=5,ncol=7)
for(i in 1:5){
  profiles_d[i,]=health_d[i,]/row_sum_d[i]
}

## Matriz quitando el centroide
R_d=matrix(0,nrow=5,ncol=7)
for(i in 1:5){
  R_d[i,]=profiles_d[i,]-centroide_d
}

## Encontramos la descomposición gsvd 
res_d=svd(Dr_sq_d%*%R_d%*%Dc_sq_d)

## Valores propios
values_d=res_d$d

## Matrices
U_d=res_d$u
V_d=res_d$v

N_d=solve(Dr_sq_d)%*%U_d
M_d=solve(Dc_sq_d)%*%V_d

## Coordenadas de los renglones
g=N_d[,c(1:2)]%*%diag(values_d[c(1:2)])

df_points_dual=data.frame(x=g[,1],y=g[,2])
df_labels_dual=data.frame(x=g[,1],y=g[,2])

ggplot(data=df_points_dual,aes(x=x,y=y))+
  geom_text(show.legend=F,col="red",size=3,
            aes(label = c("Muy Bueno","Bueno","Regular",
                          "Malo","Muy Malo")))+
  theme_minimal()+
  labs(x="",y="")+
  ylim(c(-.2,.2))

# Se grafican las dos
ggplot(data=df_labels,aes(x=x,y=y))+
  geom_text(show.legend = F,size=3,aes(label = c("16-24","25-34","35-44",
                                                 "45-54","55-64","65-74",
                                                 "75+")))+
  geom_text(show.legend=F,data=df_labels_dual,col='red',size=3,
            aes(label = c("Muy Bueno","Bueno","Regular","Malo","Muy Malo")))+
  theme_minimal()+
  labs(x="",y="")+
  ylim(c(-.2,.2))

# Existe un posible error al momento de graficar, por la no 
# unicidad de los valores singulares que están definidos
# salvo un signo. 

df_labels_dual=data.frame(x=-g[,1],y=g[,2])

# Se grafican las dos
ggplot(data=df_labels,aes(x=x,y=y))+
  geom_text(show.legend = F,size=3,aes(label = c("16-24","25-34","35-44",
                                                 "45-54","55-64","65-74",
                                                 "75+")))+
  geom_text(show.legend=F,data=df_labels_dual,col='red',size=3,
            aes(label = c("Muy Bueno","Bueno","Regular","Malo","Muy Malo")))+
  theme_minimal()+
  labs(x="",y="")+
  ylim(c(-.2,.2))
##############################################################################

##############################################################################
# Datos salud (método svd: algoritmo 2)

## Sumas por renglones y columnas
row_sum=apply(health,1,sum)
col_sum=apply(health,2,sum)

## Masas y centroide
row_masses=row_sum/n
centroide=col_sum/n

## Matriz diagonales de masas y centroide
Dr=diag(row_masses)
Dc=diag(centroide)

## Raiz cuadrada de matrices diagonales
Dr_sq=solve(sqrtm(Dr))
Dc_sq=solve(sqrtm(Dc))

## Transformamos a matrices las masas y centroides
row_m=as.matrix(row_masses)
centroide_m=as.matrix(centroide)

## Matriz de correspondencias
P=health/n

## Matriz estandarizada y centrada
A=Dr_sq%*%(P-(row_m%*%t(centroide_m)))%*%Dc_sq

## svd de A
res=svd(A)

## Valores singulares y vectores
values=res$d
U=res$u
V=res$v

## Coordenadas estándar
X=Dr_sq%*%U
Y=Dc_sq%*%V

## Coordenadas principales
f=X%*%diag(res$d)
g=Y%*%diag(res$d)


## Graficamos 
df=data.frame(x=f[,1],y=f[,2])
df_dual=data.frame(x=g[,1],y=g[,2])


ggplot(df,aes(x=x,y=y))+
  geom_text(size=3,show.legend=F,
            aes(label = c("16-24","25-34","35-44",
                          "45-54","55-64","65-74",
                          "75+")))+
  geom_text(size=3,show.legend=F,data=df_dual,col='red',
            aes(label = c("Muy Bueno","Bueno","Regular",
                          "Malo","Muy Malo")))+
  theme_minimal()+
  labs(x='',y='')+
  ylim(c(-.2,.2))+
  xlim(c(-.45,.77))
##############################################################################

##############################################################################
# Datos fumadores (método svd: algoritmo 2 visto en clase)
smoking=matrix(c(4,2,3,2,
                 4,3,7,4,
                 25,10,12,4,
                 18,24,33,13,
                 10,6,7,2),nrow=5,ncol=4,byrow=T)

## Total de observaciones
n=sum(smoking)

## Sumas por renglones y columnas
row_sum=apply(smoking,1,sum)
col_sum=apply(smoking,2,sum)

## Masas y centroide
row_masses=row_sum/n
centroide=col_sum/n

## Matriz diagonales de masas y centroide
Dr=diag(row_masses)
Dc=diag(centroide)

## Raiz cuadrada de matrices diagonales
Dr_sq=solve(sqrtm(Dr))
Dc_sq=solve(sqrtm(Dc))

## Transformamos a matrices las masas y centroides
row_m=as.matrix(row_masses)
centroide_m=as.matrix(centroide)

## Matriz de correspondencias
P=smoking/n

## Matriz estandarizada y centrada
A=Dr_sq%*%(P-(row_m%*%t(centroide_m)))%*%Dc_sq

## svd de A
res=svd(A)

## Valores singulares y vectores
values=res$d
U=res$u
V=res$v

## Coordenadas estándar
X=Dr_sq%*%U
Y=Dc_sq%*%V

## Coordenadas principales
f=X%*%diag(res$d)
g=Y%*%diag(res$d)


## Graficamos 
df=data.frame(x=-f[,1],y=f[,2])
df_dual=data.frame(x=-g[,1],y=g[,2])


p=ggplot(df,aes(x=x,y=y))+
  geom_point(size=1,alpha=0)+
  geom_text(size=3.5,show.legend=F,data=df,
            aes(
            label = c("Sr Manager","Jr Manager",
                          "Sr Employees","Jr Employees","Secretaries")))+
  theme_minimal()+
  xlim(-.3,.45)+
  labs(x="",y="")
p

p+geom_text(size=3.5,show.legend=F,data=df_dual,
            aes(colour=c("2","2","2","2"),
                label = c("None","Light","Medium","Heavy")))+
  theme_minimal()+
  labs(x="",y="")


## Añadimos un renglón
r_new=c(.42,.29,.2,.09)
f_new=numeric(2)
f_new[1]=sum(r_new*(-1*g[,1]))/values[1]
f_new[2]=sum(r_new*g[,2])/values[2]

## Graficamos
df=data.frame(x=c(-f[,1],f_new[1]),y=c(f[,2],f_new[2]))
df_dual=data.frame(x=-g[,1],y=g[,2])


p=ggplot(df,aes(x=x,y=y))+
  geom_point(size=1,alpha=0)+
  geom_text(size=3.5,show.legend=F,data=df,
            aes(colour=c("1","1","1","1","1","3"),
                label = c("Sr Manager","Jr Manager",
                          "Sr Employees","Jr Employees",
                          "Secretaries","Average")))+
  theme_minimal()+
  xlim(-.3,.45)+
  labs(x="",y="")


p+geom_text(size=3.5,show.legend=F,data=df_dual,
            aes(label = c("None","Light","Medium","Heavy")))+
  theme_minimal()+
  labs(x="",y="")


## Añadimos dos columnas
not_drinking=c(0,1,5,10,7)
drinking=c(11,17,46,78,18)

not_drinking=not_drinking/sum(not_drinking)
drinking=drinking/sum(drinking)

c1_new=numeric(2)
c1_new[1]=sum(not_drinking*-f[,1])/values[1]
c1_new[2]=sum(not_drinking*f[,2])/values[2]

c2_new=numeric(2)
c2_new[1]=sum(drinking*-f[,1])/values[1]
c2_new[2]=sum(drinking*f[,2])/values[2]

## Graficamos
df=data.frame(x=c(-f[,1],f_new[1]),y=c(f[,2],f_new[2]))
df_dual=data.frame(x=c(-g[,1],c1_new[1],c2_new[1]),
                   y=c(g[,2],c1_new[2],c2_new[2]))


p=ggplot(df,aes(x=x,y=y))+
  geom_point(size=1,alpha=0)+
  geom_text(size=3.5,show.legend=F,data=df,
            aes(label = c("Sr Manager","Jr Manager",
                          "Sr Employees","Jr Employees",
                          "Secretaries","Average")))+
  theme_minimal()+
  xlim(-.3,.45)+
  labs(x="",y="")


p+geom_text(size=3.5,show.legend=F,data=df_dual,
            aes(colour=c("1","1","1","1","2","2"),
                label = c("None","Light","Medium","Heavy",
                          "Do not drink","Drink")))+
  theme_minimal()+
  labs(x="",y="")
##############################################################################

##############################################################################
# Análisis de correspondencias múltiple vinos con matriz indicadora
vinos=matrix(c(1,0,0,0,1,0,1,1,0,0,1,0,0,1,0,1,0,1,0,1,0,1,
               0,1,0,1,0,1,0,0,1,1,0,0,1,0,1,0,0,1,1,0,1,0,
               0,1,1,0,0,1,0,0,1,1,0,1,0,0,1,0,0,1,1,0,1,0,
               0,1,1,0,0,1,0,0,1,1,0,1,0,0,1,0,1,0,1,0,1,0,
               1,0,0,0,1,0,1,1,0,0,1,0,0,1,0,1,1,0,0,1,0,1,
               1,0,0,1,0,0,1,1,0,0,1,0,1,0,0,1,1,0,0,1,0,1),
             nrow=6,ncol=22,byrow=T)

## Total de observaciones
n=sum(vinos)

## Sumas por renglones y columnas
row_sum=apply(vinos,1,sum)
col_sum=apply(vinos,2,sum)

## Masas y centroide
row_masses=row_sum/n
centroide=col_sum/n

## Matriz diagonales de masas y centroide
Dr=diag(row_masses)
Dc=diag(centroide)

## Raiz cuadrada de matrices diagonales
Dr_sq=solve(sqrtm(Dr))
Dc_sq=solve(sqrtm(Dc))

## Transformamos a matrices las masas y centroides
row_m=as.matrix(row_masses)
centroide_m=as.matrix(centroide)

## Matriz de correspondencias
P=vinos/n

## Matriz estandarizada y centrada
A=Dr_sq%*%(P-(row_m%*%t(centroide_m)))%*%Dc_sq

## svd de A
res=svd(A)

## Valores singulares y vectores
values=res$d
U=res$u
V=res$v

## Coordenadas estándar
X=Dr_sq%*%U
Y=Dc_sq%*%V

## Coordenadas principales
f=X%*%diag(res$d)
g=Y%*%diag(res$d)

## Graficamos modificando las coordenadas para que sea visible
set.seed(31415)
df=data.frame(x=f[,1],y=f[,2]+rnorm(6,0,.5))
df_dual=data.frame(x=g[,1],y=g[,2]+rnorm(6,0,.5))

p=ggplot(df_dual,aes(x=x,y=y))+
  geom_text(size=3.5,show.legend=F,data=df_dual,
            aes(label = c("Fruity-Y",
                          "Fruity-N",
                          "Woody-Y",
                          "Woody-S",
                          "Woody-N",
                          "Coffee-Y",
                          "Coffee-N",
                          "Red Fruit-Y",
                          "Red Fruit-N",
                          "Roasted-Y",
                          "Roasted-N",
                          "Vanilla-Y",
                          "Vanilla-S",
                          "Vanilla-N",
                          "Woody2-Y",
                          "Woody2-N",
                          "Fruity3-Y",
                          "Fruity3-N",
                          "Butter-Y",
                          "Butter-N",
                          "Woody3-Y",
                          "Woddy3-N")))+
  theme_minimal()+
  labs(x="",y="")
p

p=p+geom_text(size=3.5,show.legend=F,data=df,col="red",
            aes(label = c("W1","W2","W3","W4","W5","W6")))+
  theme_minimal()+
  labs(x="",y="")
p

## Añadimos un renglón
r_new=c(0,1,0,1,0,.5,.5,1,0,1,0,0,1,0,.5,.5,1,0,.5,.5,0,1)
f_new=numeric(2)
f_new[1]=sum(r_new*(1*g[,1]))/values[1]
f_new[2]=sum(r_new*g[,2])/values[2]

df_new=data.frame(x=f_new[1],y=f_new[2]-4)

## Graficamos
p+geom_text(size=3.5,data=df_new,aes(label="W7"),col="blue")

# Inercias sin modificar
inercias=res$d^2

#Inercia total
inercia = sum(inercias)

# Contribución a la inercia
contribucion=inercias/inercia;contribucion

# Inercia acumulada
cumsum(contribucion)
##############################################################################

##############################################################################
# Análisis de correspondencia múltiple para la encuesta
# sobre la familia y los cambios de rol (1994)
women = read_xls(here('../Datos/women5.xls'))

labels = as.matrix(women[,1])
women = women%>%
  select(-...1)

women=as.matrix(women)

## Total de observaciones
n=sum(women)

## Sumas por renglones y columnas
row_sum=apply(women,1,sum)
col_sum=apply(women,2,sum)

## Masas y centroide
row_masses=row_sum/n
centroide=col_sum/n

## Matriz diagonales de masas y centroide
Dr=diag(row_masses)
Dc=diag(centroide)

## Raiz cuadrada de matrices diagonales
Dr_sq=solve(sqrtm(Dr))
Dc_sq=solve(sqrtm(Dc))

## Transformamos a matrices las masas y centroides
row_m=as.matrix(row_masses)
centroide_m=as.matrix(centroide)

## Matriz de correspondencias
P=women/n

## Matriz estandarizada y centrada
A=Dr_sq%*%(P-(row_m%*%t(centroide_m)))%*%Dc_sq

## svd de A
res=svd(A)

## Valores singulares y vectores
values=res$d
U=res$u
V=res$v

## Coordenadas estándar
X=Dr_sq%*%U
Y=Dc_sq%*%V

## Coordenadas principales
f=X%*%diag(res$d)
g=Y%*%diag(res$d)

## Graficamos 
df=data.frame(x=f[,1],y=f[,2])
df_dual=data.frame(x=g[,1],y=g[,2])

ggplot(df_dual,aes(x=x,y=y))+
  geom_text(size=3.5,show.legend=F,data=df_dual,
            aes(label = labels))+
  theme_minimal()+
  labs(x="",y="")

# Inercias sin modificar
inercias=res$d^2;inercias

# Inercia total
inercia=sum(inercias)

# Contribución a la inercia
contribucion=inercias/inercia;contribucion

# Inercia acumulada
cumsum(contribucion)

# Inercias modificadas
inercias_c=numeric(length(inercias))

for(i in 1:length(inercias)){
  if(inercias[i]>1/4){
    inercias_c[i]=((4/3)*(sqrt(inercias[i])-(1/4)))^2
  }
}

# Greenacre inercia total
Jbar=(4/3)*(sum(inercias)-((16-4)/(4^2)));Jbar

# Contribución a la inercia
contribucion_c=inercias_c/Jbar;contribucion_c

# Inercia acumulada
cumsum(contribucion_c)
##############################################################################
