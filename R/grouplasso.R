#' Do the local adaptive grouped regularization
grouplasso <- function(data, index, weights=NULL, type="linear", maxit=1000, thresh=0.001, min.frac=0.1, nlam=20, delta=2, gamma=0.8, verbose=FALSE, momentum=1, reset=10, lambda=NULL, unpenalized=NULL) {
    X.transform <- NULL
    if (is.null(weights)) {weights = rep(1,nrow(data$x))}
    
    X = data$x
    Y = data$y
    
    #Center the Y matrix:
    wmeany = weighted.mean(as.vector(Y), weights)
    Y.centered = Y - wmeany
    
    #Find the adaptive group weights:
    adamodel = lsfit(y=as.matrix(Y.centered), x=as.matrix(X), intercept=FALSE, wt=weights)
    s2 = sum(weights * adamodel$residuals**2)/sum(weights)
    adapt = adamodel$coef  #[-1] #Dont include the intercept for the adaptive weights
    
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
    data = list(x=X, y=Y.centered)
    
    if (type == "linear") {
        Sol <- grouplassoRun(data,
                      index,
                      weights,
                      adaweights=adaweights,
                      thresh,
                      inner.iter=maxit,
                      outer.iter=maxit,
                      outer.thresh=thresh,
                      min.frac=min.frac,
                      nlam=nlam,
                      lambda=lambda,
                      gamma=gamma,
                      verbose=verbose,
                      momentum=momentum,
                      reset=reset)
        
        res = list()
        
        #Return the weighted mean to the intercept:
        beta = Sol$beta
        beta[1,] = beta[1,] + wmeany
        intercept = beta[1,]
        Sol$beta = beta
        
        res[['fitted']] = fitted = as.matrix(X) %*% beta
        res[['residuals']] = resid = sweep(fitted, 1, Y, '-')
        
        #Calculate the degrees of freedom used in estimating the coefficients.
        #See Wang and Leng, 2008 (Computational Statistics and Data Analysis (52) pg5279), for details 
        not.zip = matrix(0, nrow=0, ncol=length(lambda))
        group.df = matrix(0, nrow=0, ncol=ncol(beta))
        
        groups = unique(index)
        for (i in 1:length(groups)) {
            indx = which(index == groups[i])
            adaweight = adaweights[i]
            
            #group.df = rbind(group.df, apply(beta, 2, function(b) ifelse(!all(b[indx]==0), 1 + (length(indx)-1) * sqrt(sum(b[indx]**2)) / adaweight, 0)))
            group.df = rbind(group.df, apply(beta, 2, function(b) ifelse(!all(b[indx]==0), 1 + (length(indx)-1) * sqrt(sum(b[indx]**2)) * adaweight, 0)))
        }
        
        #res[['df']] = df = apply(group.df, 2, sum)
        #Naive df (add one for the intercept, which is estimated but does not appear in b)
        res[['df']] = df = drop(apply(beta, 2, function(b) sum(b!=0) + 1))
        res[['BIC']] = apply(resid, 2, function(x) sum(weights * x**2)) / s2 + log(sum(weights))*df
        res[['AIC']] = apply(resid, 2, function(x) sum(weights * x**2)) / s2 + 2*df
        res[['AICc']] = apply(resid, 2, function(x) sum(weights * x**2)) / s2 + 2*df + 2*df*(df+1)/(sum(weights)-df-1)
        
        Sol <- list(beta=Sol$beta, lambda=Sol$lambda, type=type, intercept=intercept, X.transform=X.transform, LS.coefs=adapt, adaweights=adaweights, weights=weights, results=res)
    }
    
    class(Sol) = "grouplasso"
    return(Sol)
}
