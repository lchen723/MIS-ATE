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
ATE_CP = NULL
GAMMA = NULL
DATA = NULL
PS = NULL
EAST = NULL
Ttrue = NULL
n = 600       # sample size    ## warnings would happen if n is "small" enough, e.g., n = 100
m = 100       # num of repitition of simulation
px = 3      # dim of parameters
pz = 3
p = px + pz

#### process of covariate

mu_X = rep(0,px)
mu_Z = rep(0,pz)

Sigma_X = matrix(0,px,px)
Sigma_Z = matrix(0,pz,pz)

sigma_X = 1
sigma_Z = 1
Sigma_e = diag(0.2,px)
rho_X = 0.5
rho_Z = 0.5

for(i in 1:px)
{
for(j in 1:pz)
{
Sigma_X[i,j] = sigma_X * rho_X^(abs(i-j))
Sigma_Z[i,j] = sigma_Z * rho_Z^(abs(i-j))
}
}

for(MIS in c(0.75,0.8,0.85,0.9,0.95)) {

for(time in 1:m) # repitition starts
{
set.seed(time)

X = mvrnorm(n, mu_X, Sigma_X, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)

Z = mvrnorm(n, mu_Z, Sigma_Z, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)




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



tauLR= NULL
tauRF= NULL
tau_true = NULL


for(time in 1:m) {


data_T = DATA[[time]][,2:(p+2)]
data_T = data.frame(cbind(data_T[,1],data_T[,2],data_T[,5]  ))

Tmis = data_T[,1]
Y = DATA[[time]][,1]
True = Ttrue[,time]

pr = 1 - pb
pcorr = (EAST[time,] - pr) / (1-pr-pr)

T = (pcorr>0.7)*1



X <- data_T[,2:3]


#logistic regression
ps_LR_naive = abs(glm(Tmis~X[,1]+X[,2])$fitted)
ps_LR = abs(glm(T~X[,1]+X[,2])$fitted)

#random forest
T <- as.factor(T)
Tmis <- as.factor(Tmis)
dataCorrect = cbind(T,X)
model<-train(y=T,x=X,method='rf');ps_rf<-abs(predict(model,dataCorrect,type="prob")[,2] - 0)
model_naive<-train(y=Tmis,x=X,method='rf');ps_rf_naive<-abs(predict(model_naive,data_T,type="prob")[,2] - 0)
ps_rf[which(ps_rf==0)]=0.001
ps_rf[which(ps_rf==1)]=0.999
ps_rf_naive[which(ps_rf_naive==1)]=0.999
ps_rf_naive[which(ps_rf_naive==0)]=0.001



Tmis = as.numeric(Tmis)-1
T = as.numeric(T)-1


## IPW
output_LR = (sum(T*Y/ps_LR) / sum(T/ps_LR)) - (sum((1-T)*Y/(1-ps_LR)) / sum((1-T)/(1-ps_LR)) )
output_rf = (sum(T*Y/ps_rf) / sum(T/ps_rf)) - (sum((1-T)*Y/(1-ps_rf)) / sum((1-T)/(1-ps_rf)) )
output_LR_naive = (sum(Tmis*Y/ps_LR_naive) / sum(Tmis/ps_LR_naive)) - (sum((1-Tmis)*Y/(1-ps_LR_naive)) / sum((1-Tmis)/(1-ps_LR_naive)) )
output_rf_naive = (sum(Tmis*Y/ps_rf_naive) / sum(Tmis/ps_rf_naive)) - (sum((1-Tmis)*Y/(1-ps_rf_naive)) / sum((1-Tmis)/(1-ps_rf_naive)) )



## AIPW
M1 <- glm(Y ~ Tmis + X[,1]+X[,2] , family = "binomial")
M2 <- glm(Y ~ T + X[,1]+X[,2] , family = "binomial")
mu1mis <- predict(M1, newdata = X[,1:2] %>% mutate(Tmis=1), type = "response")
mu0mis <- predict(M1, newdata = X[,1:2] %>% mutate(Tmis=0), type = "response")
mu1 <- predict(M2, newdata = X[,1:2] %>% mutate(T=1), type = "response")
mu0 <- predict(M2, newdata = X[,1:2] %>% mutate(T=0), type = "response")
output_LR_AIPW = mean(mu1 - mu0) + (mean(T*(Y-mu1)/ps_LR)) - (mean((1-T)*(Y-mu0)/(1-ps_LR)) )
output_rf_AIPW = mean(mu1 - mu0) + (mean(T*(Y-mu1)/ps_rf)) - (mean((1-T)*(Y-mu0)/(1-ps_rf)) )
output_LR_naive_AIPW = mean(mu1mis - mu0mis) + (mean(Tmis*(Y-mu1mis)/ps_LR_naive)) - (mean((1-Tmis)*(Y-mu0mis)/(1-ps_LR_naive)) )
output_rf_naive_AIPW = mean(mu1mis - mu0mis) + (mean(Tmis*(Y-mu1mis)/ps_rf_naive)) - (mean((1-Tmis)*(Y-mu0mis)/(1-ps_rf_naive)) )



tauLR = rbind(tauLR,c(output_LR_naive,output_LR_naive_AIPW,output_LR,output_LR_AIPW))
tauRF = rbind(tauRF,c(output_rf_naive,output_rf_naive_AIPW,output_rf,output_rf_AIPW))

tau_true = c(tau_true, (mean(Y[which(T==1)]) - mean(Y[which(T==0)])))


}


colMeans(tauLR) - mean(tau_true) -> biasLR
sqrt(colVars(tauLR)) -> seLR
biasLR^2 + seLR^2 -> mseLR
ATE_CP = rbind(ATE_CP, colMeans(t((t(tauLR)-1.96*seLR < mean(tau_true) & mean(tau_true) < t(tauLR)+1.96*seLR)*1)))
colMeans(tauRF) - mean(tau_true) -> biasRF
sqrt(colVars(tauRF)) -> seRF
biasRF^2 + seRF^2 -> mseRF
ATE_CP = rbind(ATE_CP, colMeans(t((t(tauRF)-1.96*seRF < mean(tau_true) & mean(tau_true) < t(tauRF)+1.96*seRF)*1)))


result = rbind(as.vector(t(round(cbind(biasLR,seLR,mseLR),3))),as.vector(t(round(cbind(biasRF,seRF,mseRF),3))))

ATE_EST = rbind(ATE_EST,result)

#cbind(Ttrue,T,data_T[,1])

}

ATE_EST


############################################################

ATE_EST800 = NULL
ATE_CP800 = NULL
GAMMA = NULL
DATA = NULL
PS = NULL
EAST = NULL
Ttrue = NULL
n = 800       # sample size    ## warnings would happen if n is "small" enough, e.g., n = 100
m = 1000       # num of repitition of simulation
px = 4      # dim of parameters
pz = 4
p = px + pz

#### process of covariate

mu_X = rep(0,px)
mu_Z = rep(0,pz)

Sigma_X = matrix(0,px,px)
Sigma_Z = matrix(0,pz,pz)

sigma_X = 1
sigma_Z = 1
Sigma_e = diag(0.2,px)
rho_X = 0.5
rho_Z = 0.5

for(i in 1:px)
{
for(j in 1:pz)
{
Sigma_X[i,j] = sigma_X * rho_X^(abs(i-j))
Sigma_Z[i,j] = sigma_Z * rho_Z^(abs(i-j))
}
}

for(MIS in c(0.75,0.8,0.85,0.9,0.95)) {

for(time in 1:m) # repitition starts
{
set.seed(time)

X = mvrnorm(n, mu_X, Sigma_X, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)

Z = mvrnorm(n, mu_Z, Sigma_Z, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)




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



tauLR= NULL
tauRF= NULL
tau_true = NULL


for(time in 1:m) {


data_T = DATA[[time]][,2:(p+2)]
data_T = data.frame(cbind(data_T[,1],data_T[,2],data_T[,5]  ))

Tmis = data_T[,1]
Y = DATA[[time]][,1]
True = Ttrue[,time]

pr = 1 - pb
pcorr = (EAST[time,] - pr) / (1-pr-pr)

T = (pcorr>0.7)*1



X <- data_T[,2:3]


#logistic regression
ps_LR_naive = abs(glm(Tmis~X[,1]+X[,2])$fitted)
ps_LR = abs(glm(T~X[,1]+X[,2])$fitted)

#random forest
T <- as.factor(T)
Tmis <- as.factor(Tmis)
dataCorrect = cbind(T,X)
model<-train(y=T,x=X,method='rf');ps_rf<-abs(predict(model,dataCorrect,type="prob")[,2] - 0)
model_naive<-train(y=Tmis,x=X,method='rf');ps_rf_naive<-abs(predict(model_naive,data_T,type="prob")[,2] - 0)
ps_rf[which(ps_rf==0)]=0.001
ps_rf[which(ps_rf==1)]=0.999
ps_rf_naive[which(ps_rf_naive==1)]=0.999
ps_rf_naive[which(ps_rf_naive==0)]=0.001



Tmis = as.numeric(Tmis)-1
T = as.numeric(T)-1


## IPW
output_LR = (sum(T*Y/ps_LR) / sum(T/ps_LR)) - (sum((1-T)*Y/(1-ps_LR)) / sum((1-T)/(1-ps_LR)) )
output_rf = (sum(T*Y/ps_rf) / sum(T/ps_rf)) - (sum((1-T)*Y/(1-ps_rf)) / sum((1-T)/(1-ps_rf)) )
output_LR_naive = (sum(Tmis*Y/ps_LR_naive) / sum(Tmis/ps_LR_naive)) - (sum((1-Tmis)*Y/(1-ps_LR_naive)) / sum((1-Tmis)/(1-ps_LR_naive)) )
output_rf_naive = (sum(Tmis*Y/ps_rf_naive) / sum(Tmis/ps_rf_naive)) - (sum((1-Tmis)*Y/(1-ps_rf_naive)) / sum((1-Tmis)/(1-ps_rf_naive)) )



## AIPW
M1 <- glm(Y ~ Tmis + X[,1]+X[,2] , family = "binomial")
M2 <- glm(Y ~ T + X[,1]+X[,2] , family = "binomial")
mu1mis <- predict(M1, newdata = X[,1:2] %>% mutate(Tmis=1), type = "response")
mu0mis <- predict(M1, newdata = X[,1:2] %>% mutate(Tmis=0), type = "response")
mu1 <- predict(M2, newdata = X[,1:2] %>% mutate(T=1), type = "response")
mu0 <- predict(M2, newdata = X[,1:2] %>% mutate(T=0), type = "response")
output_LR_AIPW = mean(mu1 - mu0) + (mean(T*(Y-mu1)/ps_LR)) - (mean((1-T)*(Y-mu0)/(1-ps_LR)) )
output_rf_AIPW = mean(mu1 - mu0) + (mean(T*(Y-mu1)/ps_rf)) - (mean((1-T)*(Y-mu0)/(1-ps_rf)) )
output_LR_naive_AIPW = mean(mu1mis - mu0mis) + (mean(Tmis*(Y-mu1mis)/ps_LR_naive)) - (mean((1-Tmis)*(Y-mu0mis)/(1-ps_LR_naive)) )
output_rf_naive_AIPW = mean(mu1mis - mu0mis) + (mean(Tmis*(Y-mu1mis)/ps_rf_naive)) - (mean((1-Tmis)*(Y-mu0mis)/(1-ps_rf_naive)) )



tauLR = rbind(tauLR,c(output_LR_naive,output_LR_naive_AIPW,output_LR,output_LR_AIPW))
tauRF = rbind(tauRF,c(output_rf_naive,output_rf_naive_AIPW,output_rf,output_rf_AIPW))

tau_true = c(tau_true, (mean(Y[which(T==1)]) - mean(Y[which(T==0)])))


}


colMeans(tauLR) - mean(tau_true) -> biasLR
sqrt(colVars(tauLR)) -> seLR
biasLR^2 + seLR^2 -> mseLR
ATE_CP800 = rbind(ATE_CP800, colMeans(t((t(tauLR)-1.96*seLR < mean(tau_true) & mean(tau_true) < t(tauLR)+1.96*seLR)*1)))
colMeans(tauRF) - mean(tau_true) -> biasRF
sqrt(colVars(tauRF)) -> seRF
biasRF^2 + seRF^2 -> mseRF
ATE_CP800 = rbind(ATE_CP800, colMeans(t((t(tauRF)-1.96*seRF < mean(tau_true) & mean(tau_true) < t(tauRF)+1.96*seRF)*1)))


result = rbind(as.vector(t(round(cbind(biasLR,seLR,mseLR),3))),as.vector(t(round(cbind(biasRF,seRF,mseRF),3))))

ATE_EST800 = rbind(ATE_EST800,result)

#cbind(Ttrue,T,data_T[,1])

}

ATE_EST800



#############################



Tab3 = rbind(ATE_EST,ATE_EST800)


#write.table(Tab3,"D://Tab3.csv",sep=",")


