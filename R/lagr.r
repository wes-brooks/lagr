#' Fit a model via local, adaptive grouped regularization
#'
#' \code{lagr} fits a LAGR model This method fits a local model at each location indicated by \code{fit.loc}.
#'
#' @param formula symbolic representation of the model
#' @param data data frame containing observations of all the terms represented in the formula
#' @param weights vector of prior observation weights (due to, e.g., overdispersion). Not related to the kernel weights.
#' @param family exponential family distribution of the response
#' @param coords matrix of locations, with each row giving the location at which the corresponding row of data was observed
#' @param fit.loc matrix of locations where the local models should be fitted
#' @param longlat \code{TRUE} indicates that the coordinates are specified in longitude/latitude, \code{FALSE} indicates Cartesian coordinates. Default is \code{FALSE}.
#' @param kernel kernel function for generating the local observation weights
#' @param bw bandwidth parameter
#' @param bw.type type of bandwidth - options are \code{dist} for distance (the default), \code{knn} for nearest neighbors (bandwidth a proportion of \code{n}), and \code{nen} for nearest effective neighbors (bandwidth a proportion of the sum of squared residuals from a global model)
#' @param tol.loc tolerance for the tuning of an adaptive bandwidth (e.g. \code{knn} or \code{nen})
#' @param varselect.method criterion to minimize in the regularization step of fitting local models - options are \code{AIC}, \code{AICc}, \code{BIC}, \code{GCV}
#' @param tuning logical indicating whether this model will be used to tune the bandwidth, in which case only the tuning criteria are returned
#' @param D pre-specified matrix of distances between locations
#' @param verbose print detailed information about our progress?
#' 
#' @return list containing the local models.
#' 
#' @export
lagr <- function(formula, data, family=gaussian(), weights=NULL, coords, fit.loc=NULL, tuning=FALSE, predict=FALSE, simulation=FALSE, oracle=NULL, kernel, bw=NULL, varselect.method=c('AIC','BIC','AICc', 'wAIC'), verbose=FALSE, longlat=FALSE, tol.loc=NULL, bw.type=c('dist','knn','nen'), D=NULL, lambda.min.ratio=0.001, n.lambda=50, lagr.convergence.tol=0.001, lagr.max.iter=20, jacknife=FALSE, bootstrap.index=NULL, na.action=na.omit, contrasts=NULL) {
    result = list()
    class(result) = "lagr"
    
    cl <- match.call()
    formula = eval.parent(substitute_q(formula, sys.frame(sys.parent())))
    mf = eval(lagr.parse.model.frame(formula, data, family, weights, coords, fit.loc, longlat, na.action, contrasts))

    y = mf$y
    x = mf$x
    w = mf$w
    mt = mf$mt
    coords = mf$coords
    dist = mf$dist
    max.dist = mf$max.dist
    min.dist = mf$min.dist
    family = mf$family
    
    #Set some variables that determine how we fit the model
    varselect.method = match.arg(varselect.method)
    
    if (is(bw, 'lagr.bw')) {
        bw.type = bw$bw.type
        bw = bw$bw
    } else {
        bw.type = match.arg(bw.type)
    }
    
    #Fit the model:
    result[['fits']] = lagr.dispatch(
        x=x,
        y=y,
        family=family,
        prior.weights=w,
        tuning=tuning,
        predict=predict,
        simulation=simulation,
        coords=coords,
        oracle=oracle,
        fit.loc=fit.loc,
        D=dist,
        varselect.method=varselect.method,
        verbose=verbose,
        bw=bw,
        bw.type=bw.type,
        kernel=kernel,
        min.dist=min.dist,
        max.dist=max.dist,
        tol.loc=tol.loc,
        lambda.min.ratio=lambda.min.ratio,
        n.lambda=n.lambda, 
        lagr.convergence.tol=lagr.convergence.tol,
        lagr.max.iter=lagr.max.iter,
        jacknife = jacknife,
        bootstrap.index=bootstrap.index
    )
    
    coefs = as.data.frame(t(sapply(result[['fits']], function(x) x[['model']][['results']][['big.avg']])))
    is.zero = as.data.frame(t(sapply(result[['fits']], function(x) x[['model']][['results']][['is.zero']])))
    varnames = rownames(result[['fits']][[1]]$coef)
    colnames(coefs) = varnames
    colnames(is.zero) = varnames
    
    #Store results from model fitting:
    if (!tuning) {
        result[['data']] = data
        result[['response']] = as.character(formula[[2]])
        result[['family']] = family
        result[['weights']] = w
        result[['coords']] = coords
        result[['fit.locs']] = fit.loc
        result[['longlat']] = longlat
        result[['kernel']] = kernel
        result[['bw']] = bw
        result[['bw.type']] = bw.type
        result[['varselect.method']] = varselect.method
        result[['dim']] = mf$dim
        result[['coefs']] = coefs
        result[['is.zero']] = is.zero
        
        result[['na.action']] <- attr(mf, "na.action")
        result[['contrasts']] <- attr(x, "contrasts")
        result[['xlevels']] <- .getXlevels(mt, mf)
        result[['call']] <- cl
        result[['terms']] <- mt
    }
    
    return(result)
}
