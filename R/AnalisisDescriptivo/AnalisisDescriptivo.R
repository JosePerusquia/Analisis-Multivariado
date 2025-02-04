############################################################################
# Análisis Descriptivo de Datos Multivariados
# Autor: José A. Perusquía Cortés
############################################################################

############################################################################
# Librerias                                                
library(TeachingDemos)   # Para caras de Chernoff
library(pracma)          # Para las curvas de Andrews
library(GGally)          # Para scatter plot usando ggplot2
library(corrplot)        # Graficas de correlacion
library(here)
############################################################################

############################################################################
# Función para generar curvas de Andrews en ggplot
source(here("andrewsCurves.R"))
############################################################################

# El comando str permite identificar de forma rápida el número de obs. y el
# número de variables así como su tipo.
str(iris)

############################################################################
# Diagrama de dispersión

# En R podemos generar un diagrama de dispersión a partir de la función plot
plot(iris[,-5])

# Para bases de datos etiquetadas se puede hacer uso de la función pairs
pairs(iris[,-5],col=iris$Species)

# Una mejor alternativa es a través de la paquetería GGally la cual permite:
# Generar un diagrama como el plot
ggpairs(iris)

# Generar un diagrama como el de pairs
ggpairs(iris, aes(color = Species, alpha = 0.5),
        upper = list(continuous = "points"))

# Generar un diagrama el cual contenga además la correlación
ggpairs(iris, aes(color = Species, alpha = 0.5),
        upper = list(continuous = wrap("cor", size = 2.5)))
############################################################################

############################################################################
# Diagramas de correlación

# La función corrplot nos permite de forma gráfica entender la correlación
# por pares de las variables. Para datos no etiquetados usamos
corrplot(cor(iris[,-5]),method="ellipse",tl.pos='l')

# Si tenemos acceso a las etiquetas podemos generar las correlaciones por 
# grupo, seleccionando los renglones apropiadados.
corrplot(cor(iris[which(iris$Species=="setosa"),-5]),method="ellipse",
         tl.pos='n')
corrplot(cor(iris[which(iris$Species=="versicolor"),-5]),method="ellipse",
         tl.pos='n')
corrplot(cor(iris[which(iris$Species=="virginica"),-5]),method="ellipse",
         tl.pos='n')
############################################################################

############################################################################
# Diagrama de estrellas
stars(iris[,-5],lwd = 1)
############################################################################

############################################################################
# Caras de Chernoff

# Una primera versión de las caras (no muy adecuada)
faces(iris[,-5])

# Una mejor alternativa son las generadas por faces2
faces2(iris[,-5])

# El orden de como graficamos las variables cambia las caras
faces2(iris[,c(4,3,2,1)])
############################################################################

############################################################################
# Curvas de Andrews (mi versión)

#Para datos no agrupados
andrewsCurves(as.matrix(iris[,-5]),iris[,5],is.grouped = F)

#Para datos agrupados
andrewsCurves(as.matrix(iris[,-5]),iris[,5])

# El orden importa y cambia las curvas
andrewsCurves(as.matrix(iris[,c(4,3,2,1)]),iris[,5])

# Curvas de Andrews con paqueteria pracma (coordenadas polares)
andrewsplot(as.matrix(iris[,-5]),iris[,5])

# En coordenadas cartesianas
andrewsplot(as.matrix(iris[,-5]),iris[,5],style="cart")
############################################################################


############################################################################
# Estadísticas Descriptivas

# Media 
summary(iris)
apply(iris[,-5],2,mean)
colMeans(iris[,-5])

# Media Usando la definición vista en clase
t(as.matrix(iris[,-5]))%*%rep(1,150)/150

# Media por grupo
by(iris[,-5],iris[,5],colMeans)

# Matriz de varianza y covarianza
var(iris[,-5])
cov(iris[,-5])

# Matriz de varianza y covarianza a partir de W y definición
W=as.matrix(sweep(iris[,-5],2,colMeans(iris[,-5])))
t(W)%*%W/(150-1)

# Varianza por grupo
by(iris[,-5],iris[,5],cov)

# Matriz de correlación
cor(iris[,-5])

# Correlación por grupo
by(iris[,-5],iris[,5],cor)
############################################################################

