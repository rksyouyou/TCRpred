## 

load("data/data.rda")
source("auxfuns/TCRpred.R")
source("auxfuns/utils.R")
source("auxfuns/pkg.R")

out = TCRpred(Y=data$Y,X=data$X,K = data$K,Z = data$Z,ntrain=500,maxiter=10,tol=0.01)
eva_metric(out$Y_true,out$Y_pred)

