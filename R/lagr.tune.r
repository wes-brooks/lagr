#' Estimate the bandwidth parameter for a lagr model
#'
#' \code{lagr.tune} estimates the bandwidth parameter for a LAGR model.
#'
#' This method calls \code{lagr} repeatedly via the \code{optimize} function, searching for the bandwidth that minimizes a bandwidth selection criterion. It returns the profiled value of the selection criterion at each bandwidth that is used in the evaluation.
#'
#' @param formula symbolic representation of the model
#' @param data data frame containing observations of all the terms represented in the formula
#' @param weights vector of prior observation weights (due to, e.g., overdispersion). Not related to the kernel weights.
#' @param family exponential family distribution of the response
#' @param coords matrix of locations, with each row giving the location at which the corresponding row of data was observed
#' @param fit.loc matrix of locations where the local models should be fitted
#' @param longlat \code{TRUE} indicates that the coordinates are specified in longitude/latitude, \code{FALSE} indicates Cartesian coordinates. Default is \code{FALSE}.
#' @param kernel kernel function for generating the local observation weights
#' @param bw bandwidth for the kernel
#' @param bw.type type of bandwidth - options are \code{dist} for distance (the default), \code{knn} for nearest neighbors (bandwidth a proportion of \code{n}), and \code{nen} for nearest effective neighbors (bandwidth a proportion of the sum of squared residuals from a global model)
#' @param bwselect.method criterion to minimize when tuning bandwidth - options are \code{AICc}, \code{BICg}, and \code{GCV}
#' @param range allowable range of the bandwidth
#' @param tol.bw global error tolerance for minimizing the bandwidth selection criterion
#' @param tol.loc local error tolerance for converting an adaptive bandwidth (e.g. \code{knn} or \code{nen}) to a distance
#' @param varselect.method criterion to minimize in the regularization step of fitting local models - options are \code{AIC}, \code{AICc}, \code{BIC}, \code{GCV}
#' @param resid.type type of residual to use (relevant for non-gaussian response) - options are \code{deviance} and \code{pearson}
#' @param tuning logical indicating whether this model will be used to tune the bandwidth, in which case only the tuning criteria are returned
#' @param a pre-specified matrix of distances between locations
#' @param verbose print detailed information about our progress?
#'
#' @return \code{list(bw, trace)} where \code{bw} minimizes the bandwidth selection criterion and trace is a data frame of each bandwidth that was tried during the optimization, along with the resulting degrees of freedom used inthe LAGR model and the value of the bandwidth selection criterion.
#' 
#' @export
lagr.tune = function(formula, data=list(), family, range=NULL, weights=NULL, coords, oracle=NULL, kernel=NULL, bw.type=c('dist','knn','nen'), varselect.method=c('AIC','BIC','AICc'), verbose=FALSE, longlat=FALSE, tol.loc=.Machine$double.eps^0.25, tol.bw=.Machine$double.eps^0.25, bwselect.method=c('AICc','GCV','BICg'), resid.type=c('deviance','pearson')) {
    if (is.null(longlat) || !is.logical(longlat)) 
        longlat <- FALSE
    if (missing(coords)) 
        stop("Observation coordinates have to be given")
    
    mf <- match.call(expand.dots = FALSE)
    #m <- match(c("formula", "data", "weights"), names(mf), 0)
    m <- match(c("formula", "data"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    dp.n <- length(model.extract(mf, "response"))
    #weights <- as.vector(model.extract(mf, "weights"))
    if (!is.null(weights) && !is.numeric(weights)) 
        stop("'weights' must be a numeric vector")
    if (is.null(weights)) 
        weights <- rep(as.numeric(1), dp.n)
    if (any(is.na(weights))) 
        stop("NAs in weights")
    if (any(weights < 0)) 
        stop("negative weights")
    y <- model.extract(mf, "response")
    
    bw.type = match.arg(bw.type)
    varselect.method = match.arg(varselect.method)
    bwselect.method = match.arg(bwselect.method)
    resid.type = match.arg(resid.type)
    
    if (!is.null(range)) {
        beta1 = min(range)
        beta2 = max(range)
    } else {
        if (bw.type == "dist") {
            bbox <- cbind(range(coords[, 1]), range(coords[, 2]))
            difmin <- spDistsN1(bbox, bbox[2, ], longlat)[1]
            if (any(!is.finite(difmin))) 
                difmin[which(!is.finite(difmin))] <- 0
            beta1 <- difmin/1000
            beta2 <- 10*difmin
        } else if (bw.type == 'knn') {
            beta1 <- 0
            beta2 <- 1
        } else if (bw.type == 'nen') {
            if (family=='binomial') {beta2 = sum(weights/(mean(y)*(1-mean(y))) * (y-mean(y))**2)}
            else if (family=='poisson') {beta2 = sum(weights/(mean(y)) * (y-mean(y))**2)}
            else if (family=='gaussian') {beta2 = sum(weights * (y-mean(y))**2)}
            beta1 = beta2/1000
        }
    }
    
    #Create a new environment, in which we will store the likelihood trace from bandwidth selection.
    oo = new.env()
    opt <- optimize(
        lagr.tune.bw,
        interval=c(beta1, beta2),
        tol=tol.bw,
        maximum=FALSE,
        formula=formula,
        coords=coords,
        env=oo,
        oracle=oracle,
        family=family,
        varselect.method=varselect.method,
        kernel=kernel,
        verbose=verbose, longlat=longlat,
        data=data,
        bw.type=bw.type,
        weights=weights,
        tol.loc=tol.loc,
        resid.type=resid.type,
        bwselect.method=bwselect.method
    )
    trace = oo$trace[!duplicated(oo$trace[,1]),]
    rm(oo)
    
    bdwt <- opt$minimum
    res <- bdwt
    return(list(bw=res, trace=trace, bwselect.method=bwselect.method, resid.type=resid.type))
}
