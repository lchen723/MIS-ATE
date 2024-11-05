read.table("D:\\NHEFS.csv",header=T, sep = ",")->data_original

library(caret)
library(survey)
library(twang)
library(gbm)
library(ncvreg)
library(MASS)
library(glmnet)
library(nleqslv)
library(resample)
library(dplyr)


B = 10
#pb = 0.95
DATA = NULL
PS = NULL
ATE_EST = NULL

data = cbind(data_original$wt82, data_original$qsmk, 
log(data_original$sbp - 50), log(data_original$dbp), log(data_original$cholestero),
 data_original$price82,  (data_original$ht),  data_original$age, data_original$sex, data_original$nerves, data_original$hbpmed, data_original$race )

data[complete.cases(data), ]->data


p = dim(data)[2]-2

n = dim(data)[1]



##########################################
## Bootstrap

for(MIS in c(0,0.75,0.8,0.85,0.9,0.95)) {

tauLR= NULL
tauRF= NULL

for(m in 1:B) {

set.seed(m)

boot = sample(n,n,replace=TRUE)

DATA[[m]] = data[boot,]



data_T = DATA[[m]][,2:(p+2)]

data_T = data.frame(data_T)

Y = DATA[[m]][,1]



Tmis <- data_T[,1]
X <- data_T[,2:(p+1)]


ps_LR_naive = abs(glm(Tmis~X[,1]+X[,2])$fitted)

Tmis <- as.factor(Tmis)
model_naive<-train(y=Tmis,x=X,method='rf');ps_rf_naive<-abs(predict(model_naive,data_T,type="prob")[,2] - 0)


pr = 1 - MIS
pcorr = abs((ps_rf_naive - pr) / (1-pr-pr))
pcorr[which(pcorr>1)]=0.999
T = rbinom(n,1,pcorr)
T <- as.factor(T)

pcorrLR = abs((ps_LR_naive - pr) / (1-pr-pr))
TLR = rbinom(n,1,pcorrLR)


#Logistic reg
ps_LR = abs(glm(TLR~X[,1]+X[,2])$fitted)

#random forest
model<-train(y=T,x=X,method='rf');ps_rf<-abs(predict(model,data_T,type="prob")[,2] - 0.001)

T = as.numeric(T)-1
Tmis = as.numeric(Tmis)-1

output_rf = (sum(T*Y/ps_rf) / sum(T/ps_rf)) - (sum((1-T)*Y/(1-ps_rf)) / sum((1-T)/(1-ps_rf)) )
output_LR = (sum(TLR*Y/ps_LR) / sum(TLR/ps_LR)) - (sum((1-TLR)*Y/(1-ps_LR)) / sum((1-TLR)/(1-ps_LR)) )




M1 <- glm(Y ~ TLR + X[,1]+X[,2] , family = "gaussian")
M2 <- glm(Y ~ T + X[,1]+X[,2] , family = "gaussian")
mu1LR <- predict(M1, newdata = X[,1:2] %>% mutate(TLR=1), type = "response")
mu0LR <- predict(M1, newdata = X[,1:2] %>% mutate(TLR=0), type = "response")
mu1 <- predict(M2, newdata = X[,1:2] %>% mutate(T=1), type = "response")
mu0 <- predict(M2, newdata = X[,1:2] %>% mutate(T=0), type = "response")
output_LR_AIPW = mean(mu1LR - mu0LR) + (mean(TLR*(Y-mu1LR)/ps_LR)) - (mean((1-TLR)*(Y-mu0LR)/(1-ps_LR)) )
output_rf_AIPW = mean(mu1 - mu0) + (mean(T*(Y-mu1)/ps_rf)) - (mean((1-T)*(Y-mu0)/(1-ps_rf)) )


tauLR = rbind(tauLR,c(output_LR,output_LR_AIPW))
tauRF = rbind(tauRF,c(output_rf,output_rf_AIPW))

}


colMeans(tauLR) -> estLR
sqrt(colVars(tauLR)) -> seLR
zLR = estLR/seLR
2*pnorm(-abs(zLR)) ->pvalueLR
colMeans(tauRF)  -> estRF
sqrt(colVars(tauRF)) -> seRF
zRF = estRF/seRF
2*pnorm(-abs(zRF)) ->pvalueRF


result = rbind(as.vector(t(round(cbind(estLR,seLR,pvalueLR),3))),as.vector(t(round(cbind(estRF,seRF,pvalueRF),3))))

ATE_EST = rbind(ATE_EST,result)



}

ATE_EST
colnames(ATE_EST) = c("EST", "SE", "p-value", "EST", "SE", "p-value")






