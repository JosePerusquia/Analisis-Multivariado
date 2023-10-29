##############################################################################
# Análisis de correspondencias
# Autor: José A. Perusquía Cortés
##############################################################################

##############################################################################
# Librerías
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(expm)
library(ca)
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
n=sum(health)

## Suma por renglones y columnas
row_sum=apply(health,1,sum)
col_sum=apply(health,2,sum)

## Masas y centroide
row_masses=row_sum/n
centroide=col_sum/n

## Matrices diagonales de masas y centroide
Dr=diag(row_masses)
Dc=diag(1/centroide)

## Raiz cuadrada de matrices diagonales
Dr_sq=sqrtm(Dr)
Dc_sq=sqrtm(Dc)

## Perfiles por renglón
profiles=matrix(0,nrow=7,ncol=5)
for(i in 1:7){
  profiles[i,]=health[i,]/row_sum[i]
}

## Matriz quitando el centroide
R=matrix(0,nrow=7,ncol=5)
for(i in 1:7){
  R[i,]=profiles[i,]-centroide
}

## Encontramos la descomposición gsvd 
res=svd(Dr_sq%*%R%*%Dc_sq)

## Valores propios
values=res$d

## Matrices
U=res$u
V=res$v

N=solve(Dr_sq)%*%U
M=solve(Dc_sq)%*%V

## Coordenadas de los renglones
f=N[,c(1:2)]%*%diag(values[c(1:2)])

df_points=data.frame(x=f[,1],y=c(0,0,0,0,0,0,0))
df_labels=data.frame(x=f[,1],y=c(.1,-.1,.1,-.1,.1,-.1,.1))

p=ggplot(data=df_points,aes(x=x,y=y))+
  geom_point(size=1)+
  geom_text(data=df_labels,show.legend=F,
            aes(label = c("16-24","25-34","35-44",
                          "45-54","55-64","65-74","75+")))+
  theme_minimal()+
  labs(x="",y="")+
  ylim(c(-.5,.5))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())
p
##############################################################################

##############################################################################
# Problema dual estado de salud
health=t(health)

## Suma por renglones y columnas
row_sum=apply(health,1,sum)
col_sum=apply(health,2,sum)

## Masas y centroide
row_masses=row_sum/n
centroide=col_sum/n

## Matrices diagonales de masas y centroide
Dr=diag(row_masses)
Dc=diag(1/centroide)

## Raiz cuadrada de matrices diagonales
Dr_sq=sqrtm(Dr)
Dc_sq=sqrtm(Dc)

## Perfiles por renglón
profiles=matrix(0,nrow=5,ncol=7)
for(i in 1:5){
  profiles[i,]=health[i,]/row_sum[i]
}

## Matriz quitando el centroide
R=matrix(0,nrow=5,ncol=7)
for(i in 1:5){
  R[i,]=profiles[i,]-centroide
}

## Encontramos la descomposición gsvd 
res=svd(Dr_sq%*%R%*%Dc_sq)

## Valores propios
values=res$d

## Matrices
U=res$u
V=res$v

N=solve(Dr_sq)%*%U
M=solve(Dc_sq)%*%V

## Coordenadas de los renglones
g=N[,c(1:2)]%*%diag(values[c(1:2)])


df_points_dual=data.frame(x=g[,1],y=c(0,0,0,0,0))
df_labels_dual=data.frame(x=g[,1],y=c(.15,-.15,.15,-.15,.15))

p=ggplot(data=df_points_dual,aes(x=x,y=y))+
  geom_point(size=1,col="red")+
  geom_text(data=df_labels_dual,show.legend=F,col="red",
            aes(label = c("Muy Bueno","Bueno","Regular","Malo","Muy Malo")))+
  theme_minimal()+
  labs(x="",y="")+
  ylim(c(-.5,.5))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())
p


p+geom_point(size=1,data=df_points)+
  geom_text(show.legend=F,data=df_labels,
            aes(label = c("16-24","25-34","35-44",
                          "45-54","55-64","65-74","75+")))+
  theme_minimal()+
  labs(x="",y="")+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())
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
# Análisis de correspondencias múltiple vinos
vinos=matrix(c(1,0,0,0,1,0,1,1,0,0,1,0,0,1,0,1,0,1,0,1,0,1,
               0,1,0,1,0,1,0,0,1,1,0,0,1,0,1,0,0,1,1,0,1,0,
               0,1,1,0,0,1,0,0,1,1,0,1,0,0,1,0,0,1,1,0,1,0,
               0,1,1,0,0,1,0,0,1,1,0,1,0,0,1,0,1,0,1,0,1,0,
               1,0,0,0,1,0,1,1,0,0,1,0,0,1,0,1,1,0,0,1,0,1,
               1,0,0,1,0,0,1,1,0,0,1,0,1,0,0,1,1,0,0,1,0,1),
             nrow=6,
             ncol=22,byrow=T)

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

## Graficamos 
df=data.frame(x=f[,1],y=f[,2])
set.seed(5)
df_dual=data.frame(x=g[,1]+rep(c(.2,-.2),11),y=g[,2]+rnorm(22,0,.15))


p=ggplot(df_dual,aes(x=x,y=y))+
  geom_point(size=1,alpha=0)+
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
p+geom_text(size=3.5,data=df_new,aes(label="W7"),col="blue")+
  ylim(-1,2.5)

# Inercias sin modificar
inercias=res$d^2

# Contribución a la inercia
contribucion=inercias/(sum(inercias))

# Inercia acumulada
cumsum(contribucion)

# Inercias modificadas
inercias_c=numeric(length(inercias))

for(i in 1:length(inercias)){
  if(inercias[i]>1/10){
    inercias_c[i]=((10/9)*(inercias[i]-(1/10)))^2
  }
}

# Greenacre inercia total
Jbar=(10/9)*(sum(inercias^2)-((22-10)/(10^2)))

# Contribución a la inercia
contribucion_c=inercias_c/Jbar

# Inercia acumulada
cumsum(contribucion_c)

##############################################################################

