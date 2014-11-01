#' Evaluate the bandwidth selection criterion for a given bandwidth
#' 
#' @param formula symbolic representation of the model
#' @param data data frame containing observations of all the terms represented in the formula
#' @param weights vector of prior observation weights (due to, e.g., overdispersion). Not related to the kernel weights.
#' @param family exponential family distribution of the response
#' @param bw bandwidth for the kernel
#' @param kernel kernel function for generating the local observation weights
#' @param coords matrix of locations, with each row giving the location at which the corresponding row of data was observed
#' @param longlat \code{TRUE} indicates that the coordinates are specified in longitude/latitude, \code{FALSE} indicates Cartesian coordinates. Default is \code{FALSE}.
#' @param varselect.method criterion to minimize in the regularization step of fitting local models - options are \code{AIC}, \code{AICc}, \code{BIC}, \code{GCV}
#' @param tol.loc tolerance for the tuning of an adaptive bandwidth (e.g. \code{knn} or \code{nen})
#' @param bw.type type of bandwidth - options are \code{dist} for distance (the default), \code{knn} for nearest neighbors (bandwidth a proportion of \code{n}), and \code{nen} for nearest effective neighbors (bandwidth a proportion of the sum of squared residuals from a global model)
#' @param bwselect.method criterion to minimize when tuning bandwidth - options are \code{AICc}, \code{BICg}, and \code{GCV}
#' @param resid.type type of residual to use (relevant for non-gaussian response) - options are \code{deviance} and \code{pearson}
#' @param verbose print detailed information about our progress?
#' 
#' @return value of the \code{bwselect.method} criterion for the given bandwidth
#' 
lagr.tune.bw = function(x, y, weights, coords, dist, family, bw, kernel, env, oracle, varselect.method, tol.loc, bw.type, bwselect.method, min.dist, max.dist, resid.type, lambda.min.ratio, n.lambda, lagr.convergence.tol, lagr.max.iter, verbose) {    
    #Fit the model with the given bandwidth:
    cat(paste("starting bw:", round(bw, 3), '\n', sep=''))
    
    lagr.model = lagr.dispatch(
        x=x,
        y=y,
        coords=coords,
        fit.loc=NULL,
        D=dist,
        family=family,
        prior.weights=weights,
        tuning=TRUE,
        predict=FALSE,
        simulation=FALSE,
        oracle=oracle,
        varselect.method=varselect.method,
        verbose=verbose,
        bw=bw,
        bw.type=bw.type,
        kernel=kernel,
        min.dist=min.dist,
        max.dist=max.dist,
        tol.loc=tol.loc,
        resid.type=resid.type,
        lambda.min.ratio=lambda.min.ratio,
        n.lambda=n.lambda, 
        lagr.convergence.tol=lagr.convergence.tol,
        lagr.max.iter=lagr.max.iter
    )
    
    # Compute the model degrees of freedom
    n = nrow(x)
    df = sum(sapply(lagr.model, function(x) tail(x[['tunelist']][['df-local']], 1)))
    dispersion = sum(sapply(lagr.model, function(x) x[['tunelist']][['dispersion']]))

    if (bwselect.method != 'jacknife') {
        fitted = sapply(lagr.model, function(x) {
            crit = x[['tunelist']][['criterion']]
            crit.weights = exp(-0.5*(min(crit)-crit)**2)
            sum(x[['tunelist']][['localfit']] * crit.weights) / sum(crit.weights)
            })
        dev.resids = family$dev.resids(y, fitted, weights)
        ll = family$aic(y, n, fitted, weights, sum(dev.resids))
    
        #Compute the loss at this bandwidth
        if (bwselect.method=='AICc') {
            loss = ll + 2*df + 2*df*(df+1)/(n-df-1)
        } else if (bwselect.method=='AIC') {
            loss = ll + 2*df
        } else if (bwselect.method=='GCV') {
            loss = ll
        } else if (bwselect.method=='CV') {
        
            loss = ll + 2*df + 2*df*(df+1)/(n-df-1)
        } else if (bwselect.method=='BICg') {
            loss = ll + log(n)*df
            #trH = sum(sapply(lagr.model, function(x) {
            #    dispersion = x[['tunelist']][['dispersion']]
            #    if (family=='gaussian') { ll = min(x[['tunelist']][['ssr-loc']][[resid.type]]) / dispersion }
            #    else if (family=='binomial') { ll = min(x[['tunelist']][['ssr-loc']][[resid.type]]) }
            #    else if (family=='poisson') { ll = min(x[['tunelist']][['ssr-loc']][[resid.type]]) /  }
            #    df = x[['tunelist']][['df']]
            #    return(ll + log(x[['tunelist']][['n']]) * df / x[['tunelist']][['n']])
            #}))
            #loss = trH + sum(sapply(lagr.model, function(x) min(x[['tunelist']][['ssr-loc']][[resid.type]])))
            #"Simplistic" BIC - based on eq4.22 from the Fotheringham et al. book:
            #loss = nrow(x) * (log(mean(sapply(lagr.model, function(x) {x[['ssr.local']]}))) + 1 + log(2*pi)) + trH * log(nrow(x))/2
        }
    } else {
        fitted = sapply(1:n, function(k) sum(lagr.model[[k]][['coefs']] * cbind(1, x[k,])))
        dev.resids = family$dev.resids(y, fitted, weights)
        loss = sum(family$aic(y, n, fitted, weights, sum(dev.resids)))
    }
    
    res = mget('trace', env=env, ifnotfound=list(matrix(NA, nrow=0, ncol=3)))
    res$trace = rbind(res$trace, c(bw, loss, df))
    assign('trace', res$trace, env=env)
    
    cat(paste('Bandwidth: ', round(bw, 3), '. df: ', round(df,4), '. Loss: ', signif(loss, 5), '\n', sep=''))
    return(loss)
}
