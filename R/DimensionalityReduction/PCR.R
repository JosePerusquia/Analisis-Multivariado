###########################################################
# Principal components regression (PCR)                                       
# Author: Jose Antonio Perusquia Cortes
# Afil : Facultad de Ciencias - UNAM
# Module : Multivariate Analysis
###########################################################

###########################################################
# Required libraries
library(ggplot2)            # Version 3.5.2
library(ggthemes)           # Version 5.1.0
library(nortest)            # Version 1.0-4
library(lmtest)             # Version 0.9-40
library(MASS)               # Version 7.3-60
library(GGally)             # Version 2.2.1
library(dplyr)              # Version 1.1.4
library(here)               # Version 1.0.1
library(factoextra)         # Version 1.0.7
library(corrplot)           # Version 0.92
library(car)                # Version 3.1-2
library(psych)              # Version 2.4.3
###########################################################

###########################################################
# Example: Customer life cycle
crm=read.csv(here("R/DimensionalityReduction/Data/dataCustomers.csv"),
             header = T,row.names=1)

# Data description and correlation plot
describe(crm)
crm%>%cor()%>%corrplot(tl.pos='l',tl.col='black',tl.cex=.75)
###########################################################

###########################################################
# Multiple linear regression to explain customer
# satisfaction
crm=crm%>%scale()%>%as.data.frame()
mod1 = lm(customerSatis~.,crm)
summary(mod1)

# Residual analysis
residuals = data.frame(x=mod1$fitted.values,
                       y=mod1$residuals)

# Normality
ggplot(data=residuals,aes(x=y))+
  geom_histogram(col='black',fill='maroon4')+
  labs(x='',y='')+
  theme_minimal()

ggplot(data=residuals,aes(x=y))+
  geom_boxplot(col='black',fill='maroon4')+
  labs(x='',y='')+
  theme_minimal()

ggplot(data=residuals,aes(sample=y))+
  geom_qq(distribution = qnorm)+
  geom_qq_line(distribution = qnorm,col='red')+
  labs(x='Theoretical quantiles',y='Empirical quantiles')+
  theme_minimal()

ad.test(residuals$y)
shapiro.test(residuals$y)

# Constant variance
ggplot(data=residuals,aes(x=x,y=y))+
  geom_point()+
  labs(x='Fitted values',y='Residuals')+
  geom_hline(yintercept =0,col='red')+
  theme_minimal()

ncvTest(mod1)

# Non-correlation test
dwtest(mod1)

# Check the multicolinearity with VIF's
vifs = vif(mod1);vifs
###########################################################

###########################################################
# Principal component analysis

# There is no need to center or scale since the data has
# been already scaled before the linear regression
PCAcrm=crm%>%select(-customerSatis)%>%
  prcomp(center=F,scale.=F)
summary(PCAcrm)

# Biplot
fviz_pca_biplot(PCAcrm,title="",ggtheme = theme_minimal(),
                geom="text", labelsize = 3)

# Variance of the components
x=c(1:16)
y=PCAcrm$sdev^2;y
cumsum(y/sum(y))

exp_var=data.frame("PC"=x,"Var"=y)
ggplot(data=exp_var,aes(x=PC,y=Var))+
  geom_col(fill="navyblue",alpha=.85)+ 
  theme_light()

# Keep the first 6 components to account for 83.20% of the
# original variability
k = 6
crm2=as.data.frame(cbind(PCAcrm$x[,c(1:k)],
                         crm$customerSatis))
names(crm2)[(k+1)]='customerSatis'

# Correlation between components and correlation with
# customer satisfaction
crm2%>%cor()%>%corrplot(tl.pos='l',tl.col='black')

# Principal component regression (PCR)
mod2 = lm(customerSatis~.-1,crm2)
summary(mod2)

# Multicollinearity is removed by transforming predictors 
# into orthogonal components, at the cost of 
# interpretability
vif(mod2)

# Residuals analysis show that normality and constant
# variance are still not met
residuals_PCA = data.frame(x=mod2$fitted.values,
                           y=mod2$residuals)

# Normality
ggplot(data=residuals_PCA,aes(x=y))+
  geom_histogram(breaks=hist(mod2$residuals,plot=F)$breaks,
                 col='black',
                 fill='maroon4')+
  labs(x='',y='')+
  theme_minimal()

ggplot(data=residuals_PCA,aes(x=y))+
  geom_boxplot(col='black',fill='maroon4')+
  labs(x='',y='')+
  theme_minimal()

ggplot(data = residuals_PCA,aes(sample=y))+
  geom_qq(distribution = qnorm)+
  geom_qq_line(distribution = qnorm,col='red')+
  labs(x='',y='')+
  theme_minimal()

ad.test(residuals_PCA$y)
shapiro.test(residuals_PCA$y)

# Constant variance
ggplot(data=residuals_PCA,aes(x=x,y=y))+
  geom_point()+
  labs(x='Fitted values',y='Residuals')+
  geom_hline(yintercept =0,col='red')+
  theme_minimal()

ncvTest(mod2)

# Correlation test
dwtest(mod2)

# Compare AIC and BIC of the models
AIC(mod1,mod2)
BIC(mod1,mod2)

# Mean squared error with 10-fold cross validation shows
# OLS predicts better than PCR
set.seed(314159)
K = 10
folds = sample(rep(1:K, length.out = nrow(crm)))

cv_error_ols = numeric(K)
cv_error_pcr = numeric(K)

for(k in 1:K){
  
  train = crm[folds != k, ]
  test  = crm[folds == k, ]
  
  # OLS
  fit1 = lm(customerSatis ~ ., data = train)
  pred1 = predict(fit1, newdata = test)
  cv_error_ols[k] = mean((test$customerSatis - pred1)^2)
  
  # PCR
  pca_train = prcomp(train %>% select(-customerSatis),
                     center = TRUE, scale. = TRUE)
  
  Z_train = pca_train$x[,1:6]
  Z_test  = scale(test %>% select(-customerSatis),
                  center = pca_train$center,
                  scale  = pca_train$scale) %*% pca_train$rotation[,1:6]
  
  fit2 = lm(train$customerSatis ~ Z_train)
  pred2 = cbind(1, Z_test) %*% coef(fit2)
  
  cv_error_pcr[k] = mean((test$customerSatis - pred2)^2)
}

mean(cv_error_ols)
mean(cv_error_pcr)
###########################################################