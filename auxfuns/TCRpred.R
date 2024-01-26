TCRpred <- function(Y,X=NULL,K=NULL,Z=NULL,sid=NULL,aaSeq=NULL,abundance=NULL,k=NULL,refm=NULL,ntrain,seed = 500,maxiter,tol){

    ##
    if(is.null(K)){
        if(is.null(sid)|is.null(aaSeq)|is.null(abundance)|is.null(refm)) stop('Error! Input data is  missing!')
        K = mk_K(sid,aaSeq,abundance,refm)
    }
    if(is.null(Z)){
        if(is.null(sid)|is.null(aaSeq)|is.null(abundance)|is.null(k)) stop('Error! Input data is  missing!')
        Z = mk_Z(sid,aaSeq,abundance,k=k)
    }

    ##
    n = nrow(K)
    top = floor(ntrain/2/log(ntrain))
    
    if(is.null(X)){
        X = matrix(1,n,ncol=1)
        colnames(X) = 'Intercept'
    }

    if(all(Y%in%c(0,1))) rtype = 'bin' else rtype = 'cont'
    set.seed(seed)
    set1 = sample(1:n,ntrain) # train
    set2 = setdiff(1:n,set1) # test
    
    ## 
    Y1 = Y[set1,,drop=FALSE]
    X1 = X[set1,,drop=FALSE]
    K1 = K[set1,set1]
    Z1 = Z[set1,]
    ## 
    Y2 = Y[set2,,drop=FALSE]
    X2 = X[set2,,drop=FALSE]
    K2 = K[set2,set1]
    Z2 = Z[set2,]
    
    ##
    sfilter = feature_screening_marginal(Y1,X1,Z1,top,prop=0.05) # feature_screening
    
    Z1 = Z1[,sfilter==1]
    Z2 = Z2[,sfilter==1]
    
    if(rtype == 'bin'){
        set.seed(seed)
        out = TCRpred_bin(Y1,X1,Z1,K1,maxiter = maxiter,tol=tol)
        Y_pred = pred_fun(X2,Z2,K2,out$beta_hat,out$gamma_hat,out$alpha_hat,'bin')
    } else {
        set.seed(seed)
        out = TCRpred_cont(Y1,X1,Z1,K1,maxiter = maxiter,tol=tol)
        Y_pred = pred_fun(X2,Z2,K2,out$beta_hat,out$gamma_hat,out$alpha_hat,'cont')
    }
    
    ##
    res = list(Y_true=Y2,Y_pred=Y_pred)
    return(res)
}


TCRpred_oneround_cont <- function(Y,X,Z,K,alpha_in){

    n = length(Y)
    ## Step 1: update response variable
    Y_upd = Y - K%*%alpha_in
    ## Step 2: feature selection
    noint_X = X[,-1,drop=FALSE] # delete intercept to avoid problem in feature selection
    res = feature_selection_lasso(Y_upd,noint_X,Z,incept=TRUE)
    filter = res$filter

    if(sum(filter)==0) {
        Z_upd = Z
    } else {
        Z_upd = Z[,filter==1,drop=FALSE]
    }
    
    lfit = res$lfit
    gamma_hat = coef(lfit)[colnames(Z_upd),,drop=FALSE]
    beta_hat = coef(lfit)[1:ncol(X),]

    ## Step 3: estimate alpha
    lambda = sqrt((ncol(X)+ncol(Z_upd))/n)
    alpha_hat = solve(diag(n)*lambda + K)%*%(Y-X%*%beta_hat-Z_upd%*%gamma_hat)

    ## Step 4: evaluate performance
    Y_pred = pred_fun(X,Z_upd,K,beta_hat,gamma_hat,alpha_hat,rtype = 'cont') # with training data 
    MSE = mean((Y - Y_pred)^2)

    res = list(MSE=MSE,Y_pred=Y_pred,beta_hat = beta_hat,gamma_hat = gamma_hat,alpha_hat = alpha_hat,Z_upd = Z_upd)
    return(res)
}

TCRpred_cont <- function(Y,X,Z,K,maxiter=10, tol=1e-2){

    n = length(Y)

    ## Initialization
    alpha_in = matrix(rep(0,n),ncol=1)
    
    MSE = mean((Y-mean(Y))^2) # mean
    MSE_diff = MSE
    iter = 0
    beta_est = NULL
    feature_mat = c(ncol(Z),sum(colnames(Z)%in%trueZ))
    
    while(iter<maxiter&MSE_diff>tol){
        iter = iter + 1
        out = TCRpred_oneround_cont(Y,X,Z,K,alpha_in)
        alpha_in = out$alpha_hat
        MSE = c(MSE,out$MSE)
        MSE_diff = abs(MSE[iter+1]-MSE[iter])
        Z = out$Z_upd
        beta_est = rbind(beta_est,out$beta_hat)
    }

    names(MSE) = paste0('iter',0:iter)
    rownames(beta_est) = paste0('iter',1:iter)

    res = list(MSE = MSE,beta_hat=out$beta_hat,gamma_hat=out$gamma_hat,alpha_hat=out$alpha_hat,beta_mat= beta_est)
    return(res)
}


TCRpred_oneround_bin <- function(Y,X,Z,K,beta_in,gamma_in,alpha_in){

    n = length(Y)
    
    ## Step 1: compute working response
    delta = X%*%beta_in+Z%*%gamma_in+K%*%alpha_in
    pi_delta = exp(delta)/(1+exp(delta))
    pi_delta[pi_delta==1] = max(max(pi_delta[pi_delta<1]),0.9999)
    pi_delta[pi_delta==0] = min(min(pi_delta[pi_delta>0]),0.0001)
    Y_wok = delta + (Y-pi_delta)/pi_delta/(1-pi_delta)

    ## Step 2: feature selection
    V = sqrt(pi_delta*(1-pi_delta))
    Omega = diag(as.numeric(V))
    Y_w = Omega%*%(Y_wok-K%*%alpha_in)
    X_w = Omega%*%X 
    Z_w = Omega%*%Z
    colnames(X_w)[1] = "weighted_intercept"
    res = feature_selection_lasso(Y_w,X_w,Z_w,incept=FALSE)
    
    filter = res$filter
    if(sum(filter)==0) {
        Z_upd = Z
    } else {
        Z_upd = Z[,filter==1,drop=FALSE]
    }
    
    lfit = res$lfit
    gamma_hat = coef(lfit)[colnames(Z_upd),,drop=FALSE]
    beta_hat = coef(lfit)[colnames(X_w),]

    ## Step 3: estimate alpha
    lambda = sqrt((ncol(X)+ncol(Z_upd))/n)
    alpha_hat = solve(diag(n)*lambda + Omega%*%K)%*%(Y_w-X_w%*%beta_hat-Omega%*%Z_upd%*%gamma_hat)

    ## Step 4: evaluate performance
    Y_pred = pred_fun(X,Z_upd,K,beta_hat,gamma_hat,alpha_hat,rtype = 'bin') # with training data 
    ## cross entropy/log loss
    ce = -mean(Y*log(Y_pred) + (1-Y)*log(1-Y_pred))

    res = list(cross_entropy=ce,Y_pred=Y_pred,beta_hat = beta_hat,gamma_hat = gamma_hat,alpha_hat = alpha_hat,Z_upd = Z_upd)
    return(res)
}



TCRpred_bin <- function(Y,X,Z,K,maxiter=10, tol=1e-2){

    n = length(Y)

    ## Initialization
    Xmat = cbind(X[,-1],Z)
    pf = c(rep(0,ncol(X)-1),rep(1,ncol(Z)))
    cvfit <- cv.glmnet(Xmat, Y, alpha = 1,penalty.factor=pf,family="binomial")
    lfit <- glmnet(Xmat, Y, alpha = 1,lambda=cvfit$lambda.min,penalty.factor=pf,family="binomial")

    gamma_in = matrix(coef(lfit)[colnames(Z),],ncol=1)
    beta_in = matrix(coef(lfit)[1:ncol(X),],ncol=1)
    alpha_in = matrix(rep(0,n),ncol=1)

    ##
    cross_entropy= Inf
    cross_entropy_diff = cross_entropy
    iter = 1
    beta_est = NULL
    

    while(iter<=maxiter&cross_entropy_diff>tol){
        out = NULL
        try({out = TCRpred_oneround_bin(Y,X,Z,K,beta_in,gamma_in,alpha_in)},silent = TRUE)#
        if(is.null(out)){
            temp_iter = 0
            while(temp_iter <=10){
                temp_iter = temp_iter + 1
                try(out = TCRpred_oneround_bin(Y,X,Z,K,beta_in,gamma_in,alpha_in),silent = TRUE)#
            }
        }
        
        if(is.null(out)){
            break
        } else {
            cross_entropy = c(cross_entropy,out$cross_entropy)
            ## only allow decrease
            cross_entropy_diff = cross_entropy[iter]-cross_entropy[iter+1]
            if(cross_entropy_diff>tol){ ## only allow decrease
                iter = iter + 1
                Z = out$Z_upd
                beta_in = out$beta_hat
                gamma_in = out$gamma_hat
                alpha_in = out$alpha_hat
                beta_est = rbind(beta_est,beta_in)
            }
        }
    }

    if(iter >= 1){
        names(cross_entropy) = paste0('iter',1:iter)
        rownames(beta_est) = paste0('iter',2:iter)
    } 

    res = list(cross_entropy= cross_entropy,beta_hat=out$beta_hat,gamma_hat=out$gamma_hat,alpha_hat=out$alpha_hat,beta_mat= beta_est)
    return(res)
}


pred_fun <- function(X,Z,K,beta,gamma,alpha,rtype){
    if(is.null(beta)){
        Y_pred = rep(NA,nrow(X))
    } else {
        Z = Z[,rownames(gamma),drop=FALSE]
        mu = X%*%beta + Z%*%gamma + K%*%alpha
        if(rtype == 'bin'){
            pmu = exp(mu)/(exp(mu)+1)
            pmu[which(exp(mu) == Inf)] = 1
            pmu[pmu==1] = max(max(pmu[pmu<1]),0.9999)
            pmu[pmu==0] = min(min(pmu[pmu>0]),0.0001)
            Y_pred = pmu
        } else {
            Y_pred = mu
        }
    }
    return(Y_pred)
}



blin <- function(X,Z,beta,gamma,Y){
    mu_est = X%*%beta + Z%*%gamma
    pmu_est = exp(mu_est)/(exp(mu_est)+1)
    pmu_est[pmu_est==1] = max(pmu_est[pmu_est<1])
    pmu_est[pmu_est==0] = min(pmu_est[pmu_est>0])
    Y_lin = mu_est + (Y-pmu_est)/(pmu_est*(1-pmu_est))
    return(Y_lin)
}

feature_selection_lasso <- function(Y,X,Z,incept=TRUE){
    Xmat = cbind(X,Z) 
    pf = c(rep(0,ncol(X)),rep(1,ncol(Z)))
    cvfit <- cv.glmnet(Xmat, Y, alpha = 1,penalty.factor=pf,intercept=incept)
    lfit <- glmnet(Xmat, Y, alpha = 1,lambda=cvfit$lambda.min,penalty.factor=pf,intercept=incept)
    filter = as.numeric(coef(lfit)[colnames(Z),]!=0)
    return(list(filter=filter,lfit=lfit))
}



feature_screening_marginal <- function(Y,X,Z,top=50,prop=0.2){

    n = length(Y)
    p = ncol(Z)
    filter = rep(0,p)
    
    ## first selection: freq
    ## at least n*prop subjects contain the 2-AA seq
    filter[colSums(Z != 0)>=n*prop] = 1

    if(all(Y%in%c(0,1))) rtype = 'bin' else rtype = 'cont'
    ## second selection: significance
    pval = rep(NA,p)
    set1 = which(filter==1)
    
    for(i in set1){
        if(rtype == 'cont')
            fit = lm(Y~0+Z[,i]+X)
        else
            fit = glm(Y~0+Z[,i]+X,family=binomial(link = "logit"))
        pval[i] = summary(fit)$coef[1,4]
    }

    set2 = which(rank(pval)<= min(top,length(set1)))
    filter[-set2] = 0
    names(filter) = colnames(Z)

    return(filter)
}

