##########################################
#######  ML for PS estimation  ###########
##########################################
 
library(caret)
library(ncvreg)
library(MASS)
library(glmnet)
library(resample)
#### General Settings
ATE_EST = NULL
GAMMA = NULL
DATA = NULL
PS = NULL
EAST = NULL
Ttrue = NULL
SEL = NULL    ## PS1 + M1
n = 600       # sample size    ## n can be changed to 800
m = 1000       # num of repitition of simulation
px = 3      # dim of parameters ## px can be changed to 4
pz = 3                          ## pz can be changed to 4
p = px + pz
 
#### process of covariate
 
mu = rep(0,p)

 
Sigma = matrix(0,p,p)
rho = 0.5
 
for(i in 1:px)
 {
 for(j in 1:pz)
 {
 Sigma[i,j] = rho^(abs(i-j))
 }
 }
 
for(MIS in c(0.75,0.8,0.85,0.9,0.95)) {

for(time in 1:m) # repitition starts
{
set.seed(time)

W = mvrnorm(n, mu, Sigma, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
X = W[,1:px]
Z = W[,(px+1):p]


 
 
g_link = sin(X[,1]) + exp(Z[,1])  # the first and the (px+1)th confounders are informative
e = 1/(1+exp(-g_link))         # pass through an inv-logit function
#e =  1 - exp(-exp(g_link))     # complement log-log model form
T = (e>0.7)*1
Ttrue = cbind(Ttrue,T)

pb = MIS
Pi = matrix(c(pb,1-pb,1-pb,pb),2,2)
E = rbind(e,1-e)
East = Pi %*% E
#Tast = (East[1,]>0.7)*1
Tast = rbinom(n,1,East[1,])
 
EAST = rbind(EAST,East[1,])
 
 
PS = rbind(PS,as.vector(e))
 
#### M1
#Y = T + abs(X%*%gamma_X + Z%*%gamma_Z) + rnorm(n,0,1)
Y = T + X[,1] + Z[,1] + rnorm(n,0,1)
#######
#### M2
#pi = 1/(1+exp(-(T+abs(g_link))))
#Y = rbinom(n,1,pi)
#######
 
data = cbind(Y,Tast,X,Z)    # true data frame
#data = cbind(Y,T,W,Z)    # naive data frame
 
DATA[[time]] = data
 
 
 
 }
 
####    End of Data Generarion   ####
 
 
 
sel = rep(0,px+pz)
selError = rep(0,px+pz) 
selLR = rep(0,px+pz)
selLRError = rep(0,px+pz) 
 
for(time in 1:m) {
 
 
data_T = data.frame(DATA[[time]][,2:(p+2)])
 
 
Tmis = data_T[,1]
Y = DATA[[time]][,1]
True = Ttrue[,time]
 
pr = 1 - MIS
pcorr = (EAST[time,] - pr) / (1-pr-pr)

T = (pcorr>0.7)*1
 
 
 
 X <- data_T[,2:(p+1)]
 
 
#random forest
T <- as.factor(T)
Tmis <- as.factor(Tmis)
dataCorrect = cbind(T,X)
model <- train(y=T,x=X,method='rf')
modelError <- train(y=Tmis,x=X,method='rf') 
 
cbind(c(1:(px+pz)),as.matrix(varImp(model)$importance)) -> variable
cbind(c(1:(px+pz)),as.matrix(varImp(modelError)$importance)) -> variableError


sel[variable[which(variable[,2]>50),1]] = sel[variable[which(variable[,2]>50),1]] + 1
selError[variableError[which(variableError[,2]>50),1]] = selError[variableError[which(variableError[,2]>50),1]] + 1 


var_sel = ncvreg(X, T , family = "binomial",penalty="lasso" ,alpha = 1)
coef = rbind(var_sel$a0,as.matrix(var_sel$beta))  # extract coefficients at all values of lambda,  including the intercept
dev2 = deviance(var_sel)
reg.df2 = var_sel$df
obj2 = BIC(var_sel)
lambda = which.min(obj2)       # search max
coef = coef[ -1, lambda]

selLR[as.vector(which(coef!=0))] = selLR[as.vector(which(coef!=0))] + 1


var_sel = ncvreg(X, Tmis , family = "binomial",penalty="lasso" ,alpha = 1)
coef = rbind(var_sel$a0,as.matrix(var_sel$beta))  # extract coefficients at all values of lambda,  including the intercept
dev2 = deviance(var_sel)
reg.df2 = var_sel$df
obj2 = BIC(var_sel)
lambda = which.min(obj2)       # search max
coef = coef[ -1, lambda]

selLRError[as.vector(which(coef!=0))] = selLRError[as.vector(which(coef!=0))] + 1
 
 }
 
SEL = rbind(SEL, rbind(c(selLRError,selLR), c(selError,sel) ))
 
 }
 

#SEL/m

##########################################################################
##########################################################################

SEL2 = NULL  ## PS2 + M1
 
for(MIS in c(0.75,0.8,0.85,0.9,0.95)) {

for(time in 1:m) # repitition starts
{
set.seed(time)

W = mvrnorm(n, mu, Sigma, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
X = W[,1:px]
Z = W[,(px+1):p]


 
 
g_link = sin(X[,1]) + exp(Z[,1])
#e = 1/(1+exp(-g_link))         # pass through an inv-logit function
e =  1 - exp(-exp(g_link))     # complement log-log model form
T = (e>0.7)*1
Ttrue = cbind(Ttrue,T)

pb = MIS
Pi = matrix(c(pb,1-pb,1-pb,pb),2,2)
E = rbind(e,1-e)
East = Pi %*% E
#Tast = (East[1,]>0.7)*1
Tast = rbinom(n,1,East[1,])
 
EAST = rbind(EAST,East[1,])
 
 
PS = rbind(PS,as.vector(e))
 
#### M1
#Y = T + abs(X%*%gamma_X + Z%*%gamma_Z) + rnorm(n,0,1)
Y = T + X[,1] + Z[,1] + rnorm(n,0,1)
#######
#### M2
#pi = 1/(1+exp(-(T+abs(g_link))))
#Y = rbinom(n,1,pi)
#######
 
data = cbind(Y,Tast,X,Z)    # true data frame
#data = cbind(Y,T,W,Z)    # naive data frame
 
DATA[[time]] = data
 
 
 
 }
 
####    End of Data Generarion   ####
 
 
 
sel = rep(0,px+pz)
selError = rep(0,px+pz) 
selLR = rep(0,px+pz)
selLRError = rep(0,px+pz) 
 
for(time in 1:m) {
 
 
data_T = data.frame(DATA[[time]][,2:(p+2)])
 
 
Tmis = data_T[,1]
Y = DATA[[time]][,1]
True = Ttrue[,time]
 
pr = 1 - MIS
pcorr = (EAST[time,] - pr) / (1-pr-pr)

T = (pcorr>0.7)*1
 
 
 
 X <- data_T[,2:(p+1)]
 
 
#random forest
T <- as.factor(T)
Tmis <- as.factor(Tmis)
dataCorrect = cbind(T,X)
model <- train(y=T,x=X,method='rf')
modelError <- train(y=Tmis,x=X,method='rf') 
 
cbind(c(1:(px+pz)),as.matrix(varImp(model)$importance)) -> variable
cbind(c(1:(px+pz)),as.matrix(varImp(modelError)$importance)) -> variableError


sel[variable[which(variable[,2]>50),1]] = sel[variable[which(variable[,2]>50),1]] + 1
selError[variableError[which(variableError[,2]>50),1]] = selError[variableError[which(variableError[,2]>50),1]] + 1 


var_sel = ncvreg(X, T , family = "binomial",penalty="lasso" ,alpha = 1)
coef = rbind(var_sel$a0,as.matrix(var_sel$beta))  # extract coefficients at all values of lambda,  including the intercept
dev2 = deviance(var_sel)
reg.df2 = var_sel$df
obj2 = BIC(var_sel)
lambda = which.min(obj2)       # search max
coef = coef[ -1, lambda]

selLR[as.vector(which(coef!=0))] = selLR[as.vector(which(coef!=0))] + 1


var_sel = ncvreg(X, Tmis , family = "binomial",penalty="lasso" ,alpha = 1)
coef = rbind(var_sel$a0,as.matrix(var_sel$beta))  # extract coefficients at all values of lambda,  including the intercept
dev2 = deviance(var_sel)
reg.df2 = var_sel$df
obj2 = BIC(var_sel)
lambda = which.min(obj2)       # search max
coef = coef[ -1, lambda]

selLRError[as.vector(which(coef!=0))] = selLRError[as.vector(which(coef!=0))] + 1
 
 }
 
SEL2 = rbind(SEL2, rbind(c(selLRError,selLR), c(selError,sel) ))
 
 }
 

#SEL2/m


##########################################################################
##########################################################################

SEL3 = NULL  ## PS1 + M2
 
for(MIS in c(0.75,0.8,0.85,0.9,0.95)) {

for(time in 1:m) # repitition starts
{
set.seed(time)

W = mvrnorm(n, mu, Sigma, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
X = W[,1:px]
Z = W[,(px+1):p]


 
 
g_link = sin(X[,1]) + exp(Z[,1])
e = 1/(1+exp(-g_link))         # pass through an inv-logit function
#e =  1 - exp(-exp(g_link))     # complement log-log model form
T = (e>0.7)*1
Ttrue = cbind(Ttrue,T)

pb = MIS
Pi = matrix(c(pb,1-pb,1-pb,pb),2,2)
E = rbind(e,1-e)
East = Pi %*% E
#Tast = (East[1,]>0.7)*1
Tast = rbinom(n,1,East[1,])
 
EAST = rbind(EAST,East[1,])
 
 
PS = rbind(PS,as.vector(e))
 
#### M1
#Y = T + abs(X%*%gamma_X + Z%*%gamma_Z) + rnorm(n,0,1)
#Y = T + X[,1] + Z[,1] + rnorm(n,0,1)
#######
#### M2
pi = 1/(1+exp(-(T+abs(g_link))))
Y = rbinom(n,1,pi)
#######
 
data = cbind(Y,Tast,X,Z)    # true data frame
#data = cbind(Y,T,W,Z)    # naive data frame
 
DATA[[time]] = data
 
 
 
 }
 
####    End of Data Generarion   ####
 
 
 
sel = rep(0,px+pz)
selError = rep(0,px+pz) 
selLR = rep(0,px+pz)
selLRError = rep(0,px+pz) 
 
for(time in 1:m) {
 
 
data_T = data.frame(DATA[[time]][,2:(p+2)])
 
 
Tmis = data_T[,1]
Y = DATA[[time]][,1]
True = Ttrue[,time]
 
pr = 1 - MIS
pcorr = (EAST[time,] - pr) / (1-pr-pr)

T = (pcorr>0.7)*1
 
 
 
 X <- data_T[,2:(p+1)]
 
 
#random forest
T <- as.factor(T)
Tmis <- as.factor(Tmis)
dataCorrect = cbind(T,X)
model <- train(y=T,x=X,method='rf')
modelError <- train(y=Tmis,x=X,method='rf') 
 
cbind(c(1:(px+pz)),as.matrix(varImp(model)$importance)) -> variable
cbind(c(1:(px+pz)),as.matrix(varImp(modelError)$importance)) -> variableError


sel[variable[which(variable[,2]>50),1]] = sel[variable[which(variable[,2]>50),1]] + 1
selError[variableError[which(variableError[,2]>50),1]] = selError[variableError[which(variableError[,2]>50),1]] + 1 


var_sel = ncvreg(X, T , family = "binomial",penalty="lasso" ,alpha = 1)
coef = rbind(var_sel$a0,as.matrix(var_sel$beta))  # extract coefficients at all values of lambda,  including the intercept
dev2 = deviance(var_sel)
reg.df2 = var_sel$df
obj2 = BIC(var_sel)
lambda = which.min(obj2)       # search max
coef = coef[ -1, lambda]

selLR[as.vector(which(coef!=0))] = selLR[as.vector(which(coef!=0))] + 1


var_sel = ncvreg(X, Tmis , family = "binomial",penalty="lasso" ,alpha = 1)
coef = rbind(var_sel$a0,as.matrix(var_sel$beta))  # extract coefficients at all values of lambda,  including the intercept
dev2 = deviance(var_sel)
reg.df2 = var_sel$df
obj2 = BIC(var_sel)
lambda = which.min(obj2)       # search max
coef = coef[ -1, lambda]

selLRError[as.vector(which(coef!=0))] = selLRError[as.vector(which(coef!=0))] + 1
 
 }
 
SEL3 = rbind(SEL3, rbind(c(selLRError,selLR), c(selError,sel) ))
 
 }
 

#SEL3/m


##########################################################################
##########################################################################

SEL4 = NULL  ## PS2 + M2
 
for(MIS in c(0.75,0.8,0.85,0.9,0.95)) {

for(time in 1:m) # repitition starts
{
set.seed(time)

W = mvrnorm(n, mu, Sigma, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
X = W[,1:px]
Z = W[,(px+1):p]


 
 
g_link = sin(X[,1]) + exp(Z[,1])
#e = 1/(1+exp(-g_link))         # pass through an inv-logit function
e =  1 - exp(-exp(g_link))     # complement log-log model form
T = (e>0.7)*1
Ttrue = cbind(Ttrue,T)

pb = MIS
Pi = matrix(c(pb,1-pb,1-pb,pb),2,2)
E = rbind(e,1-e)
East = Pi %*% E
#Tast = (East[1,]>0.7)*1
Tast = rbinom(n,1,East[1,])
 
EAST = rbind(EAST,East[1,])
 
 
PS = rbind(PS,as.vector(e))
 
#### M1
#Y = T + abs(X%*%gamma_X + Z%*%gamma_Z) + rnorm(n,0,1)
#Y = T + X[,1] + Z[,1] + rnorm(n,0,1)
#######
#### M2
pi = 1/(1+exp(-(T+abs(g_link))))
Y = rbinom(n,1,pi)
#######
 
data = cbind(Y,Tast,X,Z)    # true data frame
#data = cbind(Y,T,W,Z)    # naive data frame
 
DATA[[time]] = data
 
 
 
 }
 
####    End of Data Generarion   ####
 
 
 
sel = rep(0,px+pz)
selError = rep(0,px+pz) 
selLR = rep(0,px+pz)
selLRError = rep(0,px+pz) 
 
for(time in 1:m) {
 
 
data_T = data.frame(DATA[[time]][,2:(p+2)])
 
 
Tmis = data_T[,1]
Y = DATA[[time]][,1]
True = Ttrue[,time]
 
pr = 1 - MIS
pcorr = (EAST[time,] - pr) / (1-pr-pr)

T = (pcorr>0.7)*1
 
 
 
 X <- data_T[,2:(p+1)]
 
 
#random forest
T <- as.factor(T)
Tmis <- as.factor(Tmis)
dataCorrect = cbind(T,X)
model <- train(y=T,x=X,method='rf')
modelError <- train(y=Tmis,x=X,method='rf') 
 
cbind(c(1:(px+pz)),as.matrix(varImp(model)$importance)) -> variable
cbind(c(1:(px+pz)),as.matrix(varImp(modelError)$importance)) -> variableError


sel[variable[which(variable[,2]>50),1]] = sel[variable[which(variable[,2]>50),1]] + 1
selError[variableError[which(variableError[,2]>50),1]] = selError[variableError[which(variableError[,2]>50),1]] + 1 


var_sel = ncvreg(X, T , family = "binomial",penalty="lasso" ,alpha = 1)
coef = rbind(var_sel$a0,as.matrix(var_sel$beta))  # extract coefficients at all values of lambda,  including the intercept
dev2 = deviance(var_sel)
reg.df2 = var_sel$df
obj2 = BIC(var_sel)
lambda = which.min(obj2)       # search max
coef = coef[ -1, lambda]

selLR[as.vector(which(coef!=0))] = selLR[as.vector(which(coef!=0))] + 1


var_sel = ncvreg(X, Tmis , family = "binomial",penalty="lasso" ,alpha = 1)
coef = rbind(var_sel$a0,as.matrix(var_sel$beta))  # extract coefficients at all values of lambda,  including the intercept
dev2 = deviance(var_sel)
reg.df2 = var_sel$df
obj2 = BIC(var_sel)
lambda = which.min(obj2)       # search max
coef = coef[ -1, lambda]

selLRError[as.vector(which(coef!=0))] = selLRError[as.vector(which(coef!=0))] + 1
 
 }
 
SEL4 = rbind(SEL4, rbind(c(selLRError,selLR), c(selError,sel) ))
 
 }
 

#SEL4/m


Tab1 = rbind(SEL/m,SEL2/m,SEL3/m,SEL4/m)
colnames(Tab1) = c("X1","X2","X3","X4","X5","X6","X1","X2","X3","X4","X5","X6")

#write.table(Tab1,"D://Tab1.csv",sep=",")

