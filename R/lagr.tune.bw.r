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
lagr.tune.bw = function(formula, data, weights, family, bw, kernel, coords, longlat, env, oracle, varselect.method, tol.loc, bw.type, bwselect.method, resid.type, verbose, na.action) {    
    #Fit the model with the given bandwidth:
    cat(paste("starting bw:", round(bw, 3), '\n', sep=''))
    lagr.model = lagr(formula=formula,
                      data=data,
                      family=family,
                      weights=weights,
                      tuning=TRUE,
                      coords=coords,
                      kernel=kernel,
                      oracle=oracle,
                      bw=bw,
                      varselect.method=varselect.method,
                      verbose=verbose,
                      longlat=longlat,
                      bw.type=bw.type,
                      tol.loc=tol.loc,
                      resid.type=resid.type,
                      na.action=na.action)
    
    #Compute the loss at this bandwidth
    if (bwselect.method=='AICc') {
        #trH = sum(sapply(lagr.model[['models']], function(x) {tail(x[['tunelist']][['trace.local']],1)}))
        trH = sum(sapply(lagr.model[['models']], function(x) tail(x[['tunelist']][['df-local']],1)))
        loss = nrow(data) * (log(mean(sapply(lagr.model[['models']], function(x) {x[['tunelist']][['ssr-loc']][[resid.type]]}))) + 1 + (2*(trH+1))/(nrow(data)-trH-2) + log(2*pi))
    } else if (bwselect.method=='GCV') {
        trH = sum(sapply(lagr.model[['models']], function(x) {tail(x[['tunelist']][['trace-local']],1)})) 
        loss = sum(sapply(lagr.model[['models']], function(x) {x[['tunelist']][['ssr-loc']][[resid.type]]})) / (nrow(data)-trH)**2
    } else if (bwselect.method=='BICg') {
        trH = sum(sapply(lagr.model[['models']], function(x) {
            s2 = x[['tunelist']][['s2']]
            if (family=='gaussian') { ll = min(x[['tunelist']][['ssr-loc']][[resid.type]])/s2 + log(s2) }
            else if (family=='binomial') { ll = min(x[['tunelist']][['ssr-loc']][[resid.type]]) }
            else if (family=='poisson') { ll = min(x[['tunelist']][['ssr-loc']][[resid.type]])/s2 }
            df = x[['tunelist']][['df']]
            return(ll + log(x[['tunelist']][['n']]) * df / x[['tunelist']][['n']])
        }))
        loss = trH + sum(sapply(lagr.model[['models']], function(x) {min(x[['tunelist']][['ssr-loc']][[resid.type]])}))
        #"Simplistic" BIC - based on eq4.22 from the Fotheringham et al. book:
        #loss = nrow(data) * (log(mean(sapply(lagr.model[['models']], function(x) {x[['ssr.local']]}))) + 1 + log(2*pi)) + trH * log(nrow(data))/2
    }
    
    res = mget('trace', env=env, ifnotfound=list(matrix(NA, nrow=0, ncol=3)))
    res$trace = rbind(res$trace, c(bw, loss, trH))
    assign('trace', res$trace, env=env)
    
    cat(paste('Bandwidth: ', round(bw, 3), '. df: ', round(trH,4), '. Loss: ', signif(loss, 5), '\n', sep=''))
    return(loss)
}
