################################################################################
# Factor analysis                                         
# Author: Jose Antonio Perusquia Cortes
# Afil: Facultad de Ciencias - UNAM
# Module: Multivariate Analysis
################################################################################

################################################################################
# Required libraries
library(ggplot2)            # Version 4.0.2
library(ggthemes)           # Version 5.2.0
library(psych)              # Version 2.6.1
library(here)               # Version 1.0.2
library(qgraph)             # Version 1.9.8
library(corrplot)           # Version 0.95
library(dplyr)              # Version 1.2.0
################################################################################

################################################################################
# Source function to plot the FA model
source(here('R/DimensionalityReduction/plotFA.R'))
################################################################################

################################################################################
# Example: Grades of students                                                      
grades=read.table(here("R/DimensionalityReduction/Data/calificaciones.txt"),
                  header = T)
# Correlation matrix
R=cor(grades);R

# Names for the plot
nodesNames=c('Lin','Est','Proba','Fin','Cal')
################################################################################

################################################################################
# Factor analysis with principal factor (pa) solution and without rotation

# Using one factor
efa_f1pa=fa(r=R,nfactors=1,fm="pa",rotate="none");efa_f1pa
plotFA(efa_f1pa,nodesNames = nodesNames)

# Using two factors
efa_f2pa=fa(r=R,nfactors=2,fm="pa",rotate="none");efa_f2pa
plotFA(efa_f2pa,nodesNames=nodesNames,threshold=.3)
################################################################################

################################################################################
# Factor analysis using MLE without rotation

# We need to specify the number of observations
n = 88
p = dim(R)[1]
k = 1

# Using one factor
efa_f1mle=fa(r=R,nfactors=k,n.obs=n,fm="mle",rotate="none");efa_f1mle
plotFA(efa_f1mle)

# Hypothesis test for one factor
Lambda=efa_f1mle$loadings
Psi=diag(efa_f1mle$uniquenesses)
Sigma=(Lambda%*%t(Lambda))+Psi

# Test statistic
invSigma = solve(Sigma)
logdet = determinant(invSigma %*% R, logarithm = TRUE)$modulus
minF=tr(invSigma %*% R)-logdet[1]-p;minF

n_prime=n-1-((2*p+5)/6)-((2*k)/3)
U=n_prime*minF;U

# It is the same value as the one obtained in fa() function
efa_f1mle$STATISTIC

# Compare it against a chi-squared distribution with nu degrees of freedom
# we do not reject that one factor is adequate
nu=.5*(p-k)^2-.5*(p+k)
qchisq(.95,nu)

# Using two factors
k = 2
efa_f2mle=fa(r=R,nfactors=k,n.obs=n,fm="mle",rotate="none");efa_f2mle

# Hypothesis test to see if two factors are adequate
Lambda=efa_f2mle$loadings
Psi=diag(efa_f2mle$uniquenesses)
Sigma=(Lambda%*%t(Lambda))+Psi

# Test statistic
invSigma = solve(Sigma)
logdet = determinant(invSigma %*% R, logarithm = TRUE)$modulus
minF=tr(invSigma %*% R)-logdet[1]-p;minF

n_prime=n-1-((2*p+5)/6)-((2*k)/3)
U=n_prime*minF;U

# Critical value. We can not reject the null hypothesis of two factors
nu=.5*(p-k)^2-.5*(p+k)
qchisq(.95,nu)
################################################################################

################################################################################
# Factor analysis with two factor and varimax rotation
efa_vari=fa(r=R,nfactors=2,fm="pa",rotate="varimax");efa_vari
plotFA(efa_vari,threshold = .4,nodesNames = nodesNames)

# Plot directly the rotation

# Unrotated factors
efa_df_unrot=data.frame("F1"=efa_f2pa$loadings[,1],
                        "F2"=efa_f2pa$loadings[,2],
                        "Col"=as.factor(c("O","O","C","C","C"))
                        )

p=ggplot(efa_df_unrot,aes(x=F1,y=F2,colour=Col))+
  geom_point(size=1,alpha=0)+ # used only for legend mapping
  geom_text(show.legend=F,data=efa_df_unrot,size=5,
            aes(label = c("Lin.","Est.","Prob.","Fin.","Cal.")))+
  theme_minimal()+
  guides(colour = guide_legend(title="Type", override.aes=list(alpha=1)))+
  coord_cartesian(xlim=c(0,1),ylim=c(-.5,.5))
p

# First axis
x_seq=seq(0,1,by=.1)
y1_seq=numeric(length(x_seq))
y2_seq=numeric(length(x_seq))

for(i in 1:length(x_seq)){
  y1_seq[i]=x_seq[i]*efa_vari$rot.mat[1,1]
  y2_seq[i]=x_seq[i]*efa_vari$rot.mat[2,1]
}

df_seq_x=data.frame("x1"=y1_seq,"x2"=y2_seq)
p=p+geom_line(data=df_seq_x,aes(x=x1,y=x2),col="red")
p

# Second axis
for(i in 1:length(x_seq)){
  y1_seq[i]=x_seq[i]*efa_vari$rot.mat[1,2]
  y2_seq[i]=x_seq[i]*efa_vari$rot.mat[2,2]
}

df_seq_x=data.frame("x1"=y1_seq,"x2"=y2_seq)
p=p+geom_line(data=df_seq_x,aes(x=x1,y=x2),col="red")
p

# Factors with varimax rotation
efa_df_var=data.frame("F1"=efa_vari$loadings[,1],
                      "F2"=efa_vari$loadings[,2],
                      "Col"=as.factor(c("O","O","C","C","C"))
                      )

p=ggplot(efa_df_var,aes(x=F1,y=F2,colour=Col))+
  geom_point(size=1,alpha=0)+
  geom_text(show.legend=F,data=efa_df_var,size=5,
            aes(label = c("Lin.","Est.","Prob.","Fin.","Cal.")))+
  theme_minimal()+
  guides(colour = guide_legend(title="Examen", override.aes=list(alpha=1)))
p
################################################################################

################################################################################
# Factor analysis with two factor and oblimin rotation
efa_obli=fa(r=R,nfactors=2,fm="pa",rotate="oblimin");efa_obli
plotFA(efa_obli,threshold = 0.3, oblique = T,nodesNames = nodesNames)

# Factors with oblimin rotation
efa_df_obl=data.frame("F1"=efa_obli$loadings[,1],
                      "F2"=efa_obli$loadings[,2],
                      "Col"=as.factor(c("O","O","C","C","C"))
                      )

p=ggplot(efa_df_obl,aes(x=F1,y=F2,colour=Col))+
  geom_point(size=1,alpha=0)+
  geom_text(show.legend=F,data=efa_df_obl,size=5,
            aes(label = c("Lin.","Est.","Prob.","Fin.","Cal.")))+
  theme_minimal()+
  guides(colour = guide_legend(title="Type", override.aes=list(alpha=1)))
p
################################################################################

################################################################################
# Intercorrelations of nine psychological variables
R_Holzinger = read.table(here("R/DimensionalityReduction/Data/R_Holzinger.txt"),
                         header = F)
R_Holzinger = R_Holzinger+t(R_Holzinger)
diag(R_Holzinger) = 1
nodesNamesRH = c('WM','SC','OW','MA','R','MN','G','B','H')
names(R_Holzinger) = nodesNamesRH

R_Holzinger%>%cor()%>%corrplot(tl.pos='l',tl.col='black',tl.cex=.75)

# Using one factor without rotation
efa_f1RH=fa(r=R_Holzinger,nfactors=1,fm="pa",rotate="none");efa_f1RH
plotFA(efa_f1RH,nodesNames = nodesNamesRH)

# Using two factors
efa_f2RH=fa(r=R_Holzinger,nfactors=2,fm="pa",rotate="none");efa_f2RH
plotFA(efa_f2RH,nodesNames = nodesNamesRH)

efa_f2RHV=fa(r=R_Holzinger,nfactors=2,fm="pa",rotate="varimax");efa_f2RHV
plotFA(efa_f2RHV,nodesNames = nodesNamesRH,oblique = F)
plotFA(efa_f2RHV,nodesNames = nodesNamesRH,oblique = F,threshold = .3)

efa_f2RHO=fa(r=R_Holzinger,nfactors=2,fm="pa",rotate="oblimin");efa_f2RHO
plotFA(efa_f2RHO,nodesNames = nodesNamesRH,oblique = T)
plotFA(efa_f2RHO,nodesNames = nodesNamesRH,oblique = T,threshold = .3)

# Using three factors
efa_f3RH=fa(r=R_Holzinger,nfactors=3,fm="pa",rotate="none");efa_f3RH
plotFA(efa_f3RH,nodesNames = nodesNamesRH)

efa_f3RHV=fa(r=R_Holzinger,nfactors=3,fm="pa",rotate="varimax");efa_f3RHV
plotFA(efa_f3RHV,nodesNames = nodesNamesRH,oblique = F)
plotFA(efa_f3RHV,nodesNames = nodesNamesRH,oblique = F,threshold = .3)

efa_f3RHO=fa(r=R_Holzinger,nfactors=3,fm="pa",rotate="oblimin");efa_f3RHO
plotFA(efa_f3RHO,nodesNames = nodesNamesRH,oblique = T)
plotFA(efa_f3RHO,nodesNames = nodesNamesRH,oblique = T,threshold = .3)
################################################################################



