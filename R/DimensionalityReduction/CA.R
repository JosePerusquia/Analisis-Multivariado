############################################################
# Correspondence analysis                                        
# Author: Jose Antonio Perusquia Cortes
# Afil : Facultad de Ciencias - UNAM
# Module : Multivariate analysis
############################################################

############################################################
# Required libraries 
library(ggplot2)          # Version 4.0.2
library(ggthemes)         # Version 5.2.0
library(RColorBrewer)     # Version 1.1-3
library(expm)             # Version 1.0-0
library(ca)               # Version 0.71.1
library(here)             # Version 1.0.2
library(dplyr)            # Version 1.2.1
library(readxl)           # Version 1.4.5
############################################################

############################################################
# Source functions
source('caAlgorithms.R')
############################################################

############################################################
# Example: Health status by age group
health=matrix(c(243,789,167,18,6,
                220,809,164,35,6,
                147,658,181,41,8,
                90,469,236,50,16,
                53,414,306,106,30,
                44,267,284,98,20,
                20,136,157,66,17),
              nrow=7,ncol=5,byrow=T)

# Age group
ageGroup = c("16-24","25-34","35-44","45-54","55-64",
             "65-74","75+")

# Health status
healthStatus = c("Muy Bueno","Bueno","Regular","Malo",
                 "Muy Malo")

rownames(health)=ageGroup
colnames(health)=healthStatus

# Correspondance analysis for the rows
ca_health_rows = ca_gsvd(health)

df_labels=data.frame(x=ca_health_rows[,1],
                     y=ca_health_rows[,2])

ggplot(data=df_labels,aes(x=x,y=y))+
  geom_text(show.legend=F,size=3,aes(label = ageGroup))+
  theme_minimal()+
  labs(x="",y="")+
  coord_cartesian(ylim=c(-.1,.1),xlim=c(-.45,.77))

# Correspondance analysis for the columns
ca_health_cols = ca_gsvd(t(health))

df_labels_dual=data.frame(x=ca_health_cols[,1],
                          y=ca_health_cols[,2])

ggplot(data=df_labels_dual,aes(x=x,y=y))+
  geom_text(show.legend=F,size=3,aes(label=healthStatus))+
  theme_minimal()+
  labs(x="",y="")+
  coord_cartesian(ylim=c(-.1,.1),xlim=c(-1,.6))

# Plot both the rows and the columns
ggplot(data=df_labels,aes(x=x,y=y))+
  geom_text(show.legend = F,size=3,aes(label = ageGroup))+
  geom_text(show.legend=F,data=df_labels_dual,col='red',
            size=3,aes(label = healthStatus))+
  theme_minimal()+
  labs(x="",y="")+
  coord_cartesian(ylim=c(-.1,.1),xlim=c(-1,.77))

# There's a possible issue when plotting due to the
# uniqueness of the singular values up to a sign
df_labels_dual=data.frame(x=-ca_health_cols[,1],
                          y=ca_health_cols[,2])

ggplot(data=df_labels,aes(x=x,y=y))+
  geom_text(show.legend = F,size=3,aes(label = ageGroup))+
  geom_text(show.legend=F,data=df_labels_dual,col='red',size=3,
            aes(label = healthStatus))+
  theme_minimal()+
  labs(x="",y="")+
  coord_cartesian(ylim=c(-.1,.1),xlim=c(-.45,.77))

# Second algorithm using the SVD decomposition
ca_health = ca_svd(health)

df_health_rows=data.frame(x=ca_health$f[,1],
                          y=ca_health$f[,2])
df_health_cols=data.frame(x=ca_health$g[,1],
                          y=ca_health$g[,2])

ggplot(df_health_rows,aes(x=x,y=y))+
  geom_text(size=3,show.legend=F,aes(label = ageGroup))+
  geom_text(size=3,show.legend=F,data=df_health_cols,
            col='red',aes(label = healthStatus))+
  theme_minimal()+
  labs(x='',y='')+
  coord_cartesian(ylim=c(-.1,.1),xlim=c(-.45,.8))
############################################################

############################################################
# Smoking habits SVD algorithm
smoking = matrix(c(4,2,3,2,4,3,7,4,25,10,12,4,18,24,33,
                   13,10,6,7,2),nrow=5,ncol=4,byrow=T)

# Job titles
positions = c("Sr Manager","Jr Manager","Sr Employees",
              "Jr Employees","Secretaries")

# Smoking habits
habits = c("None","Light","Medium","Heavy")

# Solution
ca_smoking = ca_svd(smoking)

df_smoking_rows=data.frame(x=ca_smoking$f[,1],
                           y=ca_smoking$f[,2])
df_smoking_cols=data.frame(x=ca_smoking$g[,1],
                           y=ca_smoking$g[,2])

ggplot(df_smoking_rows,aes(x=x,y=y))+
  geom_text(size=3,show.legend=F,aes(label = positions))+
  geom_text(size=3,show.legend=F,data=df_smoking_cols,
            col='red',aes(label = habits))+
  theme_minimal()+
  labs(x='',y='')+
  coord_cartesian(ylim=c(-.2,.15),xlim=c(-.45,.3))

# Add a new row
r_new=c(0.42,0.29,0.2,0.09)
f_new=numeric(2)
f_new[1]=sum(r_new*ca_smoking$g[,1])/ca_smoking$sv[1]
f_new[2]=sum(r_new*ca_smoking$g[,2])/ca_smoking$sv[2]

df_smoking_rows=data.frame(x=c(ca_smoking$f[,1],f_new[1]),
                           y=c(ca_smoking$f[,2],f_new[2]))

positionsMod = c(positions,'Average')
colours = c(rep('black',5),'blue')

ggplot(df_smoking_rows,aes(x=x,y=y))+
  geom_text(size=3,show.legend=F,aes(label = positionsMod),
            colour=colours)+
  geom_text(size=3,show.legend=F,data=df_smoking_cols,
            col='red',aes(label = habits))+
  theme_minimal()+
  labs(x='',y='')+
  coord_cartesian(ylim=c(-.2,.15),xlim=c(-.45,.3))

# Add two columns about drinking habits
not_drinking=c(0,1,5,10,7)
drinking=c(11,17,46,78,18)

# Normalise to work with profiles not raw counts
not_drinking=not_drinking/sum(not_drinking)
drinking=drinking/sum(drinking)

c1_new=numeric(2)
c1_new[1]=sum(not_drinking*ca_smoking$f[,1])/ca_smoking$sv[1]
c1_new[2]=sum(not_drinking*ca_smoking$f[,2])/ca_smoking$sv[2]

c2_new=numeric(2)
c2_new[1]=sum(drinking*ca_smoking$f[,1])/ca_smoking$sv[1]
c2_new[2]=sum(drinking*ca_smoking$f[,2])/ca_smoking$sv[2]

# Plot the data
df_smoking_cols=data.frame(x=c(ca_smoking$g[,1],c1_new[1],c2_new[1]),
                           y=c(ca_smoking$g[,2],c1_new[2],c2_new[2]))

habitsMod = c(habits,"Do not drink","Drink")
colours = c(rep('red',4),rep('blue',2))

ggplot(df_smoking_rows,aes(x=x,y=y))+
  geom_text(size=3,show.legend=F,aes(label = positionsMod))+
  geom_text(size=3,show.legend=F,data=df_smoking_cols,
            aes(label = habitsMod),colour=colours)+
  theme_minimal()+
  labs(x='',y='')+
  coord_cartesian(ylim=c(-.2,.35),xlim=c(-.5,.3))
############################################################

############################################################
# Example MCA: Wine and three experts 
wines=matrix(c(1,0,0,0,1,0,1,1,0,0,1,0,0,1,0,1,0,1,0,1,0,1,
               0,1,0,1,0,1,0,0,1,1,0,0,1,0,1,0,0,1,1,0,1,0,
               0,1,1,0,0,1,0,0,1,1,0,1,0,0,1,0,0,1,1,0,1,0,
               0,1,1,0,0,1,0,0,1,1,0,1,0,0,1,0,1,0,1,0,1,0,
               1,0,0,0,1,0,1,1,0,0,1,0,0,1,0,1,1,0,0,1,0,1,
               1,0,0,1,0,0,1,1,0,0,1,0,1,0,0,1,1,0,0,1,0,1),
             nrow=6,ncol=22,byrow=T)

# Id of the wines
winesId = c("W1","W2","W3","W4","W5","W6")

# Opinions by experts
opinions = c("Fruity-Y","Fruity-N","Woody-Y","Woody-S",
             "Woody-N","Coffee-Y","Coffee-N","Red Fruit-Y",
             "Red Fruit-N","Roasted-Y","Roasted-N","Vanilla-Y",
             "Vanilla-S","Vanilla-N","Woody2-Y","Woody2-N",
             "Fruity3-Y","Fruity3-N","Butter-Y","Butter-N",
             "Woody3-Y","Woddy3-N")

# Solution
ca_wines = ca_svd(wines)
df_wines_rows=data.frame(x=ca_wines$f[,1],y=ca_wines$f[,2])
df_wines_cols=data.frame(x=ca_wines$g[,1],y=ca_wines$g[,2])

ggplot(df_wines_rows,aes(x=x,y=y))+
  geom_text(size=3,show.legend=F,aes(label = winesId),
            col = 'red',position='jitter')+
  geom_text(size=3,show.legend=F,data=df_wines_cols,
            aes(label = opinions),position='jitter')+
  theme_minimal()+
  labs(x='',y='')

# Add a new row about a new wine
r_new=c(0,1,0,1,0,.5,.5,1,0,1,0,0,1,0,.5,.5,1,0,.5,.5,0,1)
f_new=numeric(2)
f_new[1]=sum(r_new*ca_wines$g[,1])/ca_wines$sv[1]
f_new[2]=sum(r_new*ca_wines$g[,2])/ca_wines$sv[2]

df_wines_rows=data.frame(x=c(ca_wines$f[,1],f_new[1]),
                         y=c(ca_wines$f[,2],f_new[2]))

winesIDmod = c(winesId,'W7')
colours = c(rep('red',6),'blue')

ggplot(df_wines_rows,aes(x=x,y=y))+
  geom_text(size=3,show.legend=F,aes(label = winesIDmod),
            col = colours,position='jitter')+
  geom_text(size=3,show.legend=F,data=df_wines_cols,
            aes(label = opinions),position='jitter')+
  theme_minimal()+
  labs(x='',y='')

# Inertia analysis
inertias=ca_wines$sv^2

# Total inertia
inertia = sum(inertias)

# Contribution to the total inertia
contribution=inertias/inertia

# Cummulative inertia
cumsum(contribution)
############################################################

############################################################
# Example: MCA for the survey about family and change of 
# roles in 1994 using Burt matrix
women = read_xls(here('R/DimensionalityReduction/Data/women5.xls'))

# Questions using Burt format
labels = as.matrix(women[,1])
women = women%>%select(-...1)
women=as.matrix(women)

# Solution 
ca_women = ca_svd(women)

df_women_rows=data.frame(x=ca_women$f[,1],y=ca_women$f[,2])
df_women_cols=data.frame(x=ca_women$g[,1],y=ca_women$g[,2])

ggplot(df_wines_rows,aes(x=x,y=y))+
  geom_text(size=3,show.legend=F,data=df_women_cols,
            aes(label = labels),position='jitter')+
  theme_minimal()+
  labs(x='',y='')

# Inertia analysis
inertias=ca_women$sv^2;inertias

# Total inertia
inertia=sum(inertias)

# Contribution to the total inertia
contribution=inertias/inertia;contribution

# Cummulate contribution
cumsum(contribution)

# Modified inertias
# Number of questions
K = 4
# Total number of levels
J = dim(women)[1]
inertias_c=numeric(length(inertias))

for(i in 1:length(inertias)){
  if(inertias[i]>1/K){
    inertias_c[i]=((K/(K-1))*(sqrt(inertias[i])-(1/K)))^2
  }
}

# Greenacre total inertia
Jbar=(K/(K-1))*(sum(inertias)-((J-K)/(K^2)))

# Contribution to the inertia
contribution_c=inertias_c/Jbar;

# Cummulative contribution
cumsum(contribution_c)
############################################################