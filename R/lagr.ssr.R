#' Calculate the sum of squared residuals for a local model
#' 
#' This function fits a local LAGR model at \code{loc}, and returns its sum of squared residuals (SSR) as a proportion of the SSR from a global model. This proportion is how the bandwidth is specified under \code{nen}.
#' 
#' @param bw kernel bandwidth (distance) to use for fitting the local model
#' @param x matrix of observed covariates
#' @param y vector of observed responses
#' @param family exponential family distribution of the response
#' @param loc location around which to center the kernel
#' @param coords matrix of locations, with each row giving the location at which the corresponding row of data was observed
#' @param dist vector of distances from central location to the observation locations
#' @param kernel kernel function for generating the local observation weights
#' @param longlat \code{TRUE} indicates that the coordinates are specified in longitude/latitude, \code{FALSE} indicates Cartesian coordinates. Default is \code{FALSE}.
#' @param bw bandwidth parameter
#' @param bw.type type of bandwidth - options are \code{dist} for distance (the default), \code{knn} for nearest neighbors (bandwidth a proportion of \code{n}), and \code{nen} for nearest effective neighbors (bandwidth a proportion of the sum of squared residuals from a global model)
#' @param tol.loc tolerance for the tuning of an adaptive bandwidth (e.g. \code{knn} or \code{nen})
#' @param varselect.method criterion to minimize in the regularization step of fitting local models - options are \code{AIC}, \code{AICc}, \code{BIC}, \code{GCV}
#' @param resid.type type of residual to use (relevant for non-gaussian response) - options are \code{deviance} and \code{pearson}
#' @param tuning logical indicating whether this model will be used to tune the bandwidth, in which case only the tuning criteria are returned
#' @param D pre-specified matrix of distances between locations
#' @param verbose print detailed information about our progress?
#' 
#' 
lagr.ssr = function(bw, x, y, family, loc, coords, dist, kernel, target, varselect.method, resid.type, prior.weights, oracle, verbose) {
    #Calculate the local weights:
    kernel.weights = drop(kernel(dist, bw))
    
    lagr.object = lagr.fit.inner(
        x=x,
        y=y,
        family=family,
        coords=coords,
        loc=loc,
        varselect.method=varselect.method,
        predict=TRUE,
        tuning=FALSE,
        simulation=FALSE,
        verbose=verbose,
        kernel.weights=kernel.weights,
        prior.weights=prior.weights,
        oracle=oracle
    )
    
    loss = lagr.object[['ssr']][[resid.type]]
    if (verbose) { cat(paste('loc:(', paste(round(loc,3), collapse=","), '), target: ', round(target,3), ', bw:', round(bw,3), ', ssr:', round(loss,3), ', miss:', round(abs(loss-target),3), '\n', sep="")) }
    
    
    return(abs(loss-target))
}
