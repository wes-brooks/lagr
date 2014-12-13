#' Do the local adaptive grouped regularization step.
#' 
#' @param data list containing \code{x}, a matrix of local covariates, and \code{y}, the corresponding vector of observed responses
#' @param index vector indicating the group membership of each column of the covariate vector
#' @param weights vector of observation weights
#' @param maxit maximum iterations to run blockwise coordinate descent
#' @param thresh iterate blockwise coordinate descent until the norm of the coefficient vector changes by less than this threshold
#' @param min.frac ratio between the smallest and largest lambdas (lasso tuning parameters)
#' @param nlam number of different lambdas (lasso tuning parameters) at which to fit the coefficients
#' @param delta exponent of the unpenalized group coefficient norm in the adaptive penalty weight
#' @param gamma 
#' @param verbose print detailed information about model fitting?
#' @param reset 
#' @param lambda vector of prespecified lasso tuning parameters - leave NULL to have the lambdas calculated automatically.
#' @param unpenalized index of any unpenalized groups
#' 
#' @return a list containing the coefficients, tuning parameters, AIC/AICc/BIC/GCV values, degrees of freedom, fitted values, and residuals
grouplasso <- function(data, index, family, weights=NULL, maxit=1000, thresh=0.001, min.frac=0.1, nlam=20, delta=2, optim.step.size=0.8, verbose=FALSE, reset=10, lambda=NULL, unpenalized=NULL) {
    X.transform <- NULL
    if (is.null(weights)) {weights = rep(1,nrow(data$x))}
    
    X = data$x
    Y = data$y
    
    #Find the adaptive group weights:
    adamodel = glm(Y~X-1, weights=weights, family=family)
    adamodel$df.residual = sum(weights)
    dispersion = summary(adamodel)$dispersion
    adapt = adamodel$coef
    
    groups = unique(index)
    n.g = length(groups)
    adaweights = rep(1, n.g)
    for (i in 1:n.g) {
        g = groups[i]
        if (g %in% unpenalized) {
            adaweights[i]=0
        } else {
            indx = which(index == g)
            adaweights[i] = 1 / sqrt(sum(adapt[indx]**2))**delta
        }    
    }
    
    #Indicate the groups whose coefficients are unpenalized:
    for (g in unpenalized) {
        indx = which(groups == g)
        adaweights[indx] = 0
    }      
    data = list(x=X, y=Y)
    
    Sol <- grouplassoRun(data,
                         index,
                         weights,
                         adaweights=adaweights,
                         family=family,
                         thresh,
                         inner.iter=maxit,
                         outer.iter=maxit,
                         outer.thresh=thresh,
                         min.frac=min.frac,
                         nlam=nlam,
                         lambda=lambda,
                         optim.step.size=optim.step.size,
                         verbose=verbose,
                         reset=reset)
    
    res = list()
    
    #Return the weighted mean to the intercept:
    beta = cbind(Sol$beta, adapt)
    intercept = drop(beta[1,])

    fitted = family$linkinv(as.matrix(X) %*% beta)
    resid = sweep(fitted, 1, Y, '-')
    dev.resid.values = apply(fitted, 2, function(mu) family$dev.resids(Y, mu, weights))
    
    #Calculate the degrees of freedom used in estimating the coefficients.
    #See Wang and Leng, 2008 (Computational Statistics and Data Analysis (52) pg5279), for details 
    not.zip = matrix(0, nrow=0, ncol=length(lambda))
    group.df = matrix(0, nrow=0, ncol=ncol(beta))
    
    groups = unique(index)
    for (i in 1:length(groups)) {
        indx = which(index == groups[i])
        adaweight = adaweights[i]
        
        group.df = rbind(group.df, apply(beta, 2, function(b) ifelse(!all(b[indx]==0), 1 + (length(indx)-1) * sqrt(sum(b[indx]**2)) * adaweight, 0)))
    }
    
    #res[['df']] = df = apply(group.df, 2, sum)
    #Naive df (add one for the intercept, which is estimated but does not appear in b)
    df = drop(apply(beta, 2, function(b) sum(b!=0) + 1))
    ll = apply(dev.resid.values, 2, sum) / dispersion
    res[['BIC']] = ll + log(sum(weights))*df
    res[['AIC']] = ll + 2*df
    res[['AICc']] = ll + 2*df + 2*df*(df+1)/(sum(weights)-df-1)

    #For each variable combination represented in the results, find the coefficients that minimise the loss:
    models = as.factor(apply(beta, 2, function(x) paste(as.integer(x!=0), collapse="")))  
    best.index = vector()
    for (l in levels(models)) {
        best.index = c(best.index, which.min(ifelse(models==l, ll, Inf)))
    }
    
    #Filter the results:
    res[['df']] = df[best.index]
    res[['BIC']] = res[['BIC']][best.index]
    res[['AIC']] = res[['AIC']][best.index]
    res[['AICc']] = res[['AICc']][best.index]
    beta = beta[,best.index]
    intercept = intercept[best.index]
    res[['fitted']] = fitted[,best.index]
    residuals = resid[,best.index]
    dev.resid.values = dev.resid.values[,best.index]
    
    #Get the AIC-optimal weights
    w.lik = sqrt(t(dev.resid.values)) %*% sqrt(dev.resid.values) / dispersion
    constraint.mat = cbind(1, diag(rep(1, length(best.index))))
    constraint.vec = c(1, rep(0,length(best.index)))
    res[['wAIC']] = -solve.QP(w.lik, -res[['df']], constraint.mat, constraint.vec, meq=1)$solution

    p.max = ncol(X)
    res[['dispersion']] = dispersion = sum(weights * dev.resid.values[,ncol(dev.resid.values)]^2) / (sum(weights)-p.max-1)
    res[['full.model.cov']] = dispersion * chol2inv(adamodel$qr$qr[1:p.max, 1:p.max])
    
    Sol <- list(beta=beta, lambda=Sol$lambda[best.index], intercept=intercept, LS.coefs=adapt, adaweights=adaweights, weights=weights, results=res, adamodel=adamodel)
    
    class(Sol) = "grouplasso"
    return(Sol)
}
