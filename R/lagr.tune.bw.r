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
#' @param verbose print detailed information about our progress?
#' 
#' @return value of the \code{bwselect.method} criterion for the given bandwidth
#' 
lagr.tune.bw = function(x, y, weights, coords, dist, family, bw, kernel, env, oracle, varselect.method, tol.loc, bw.type, bwselect.method, min.dist, max.dist, lambda.min.ratio, n.lambda, lagr.convergence.tol, lagr.max.iter, verbose) {    
    #Fit the model with the given bandwidth:
    cat(paste('Bandwidth: ', round(bw, 3), '; ', sep=""))

    # Tell lagr.dispatch whether to select bandwidth via the jacknife
    if (bwselect.method=='jacknife') {
        jacknife = TRUE
    } else {
        jacknife = FALSE
    }


    vcr.model = lagr.dispatch(
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
        tol.loc=tol.loc,,
        lambda.min.ratio=lambda.min.ratio,
        n.lambda=n.lambda, 
        lagr.convergence.tol=lagr.convergence.tol,
        lagr.max.iter=lagr.max.iter
    )
    
    res = mget('trace', env=env, ifnotfound=list(matrix(NA, nrow=0, ncol=3)))
    res$trace = as.data.frame(rbind(res$trace, c(bw, vcr.model[[bwselect.method]], vcr.model$df)))
    colnames(res$trace) = c("bw", "loss", "df")
    res$trace = res$trace[order(res$trace$bw),]
    assign('trace', res$trace, env=env)
    
    cat(paste('df: ', round(vcr.model$df,4), '; Loss: ', signif(vcr.model[[bwselect.method]], 5), '\n', sep=''))
    return(vcr.model[[bwselect.method]])
}
