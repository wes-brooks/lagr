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
#' @param resid.type type of residual to use (relevant for non-gaussian response) - options are \code{deviance} and \code{pearson}
#' @param tuning logical indicating whether this model will be used to tune the bandwidth, in which case only the tuning criteria are returned
#' @param D pre-specified matrix of distances between locations
#' @param verbose print detailed information about our progress?
#' 
#' @return list containing the local models.
#' 
#' @export
lagr <- function(formula, data, family=c('gaussian', 'binomial', 'poisson', 'Cox'), weights=NULL, coords, fit.loc=NULL, tuning=FALSE, predict=FALSE, simulation=FALSE, oracle=NULL, kernel, bw=NULL, varselect.method=c('AIC','BIC','AICc'), verbose=FALSE, longlat, tol.loc=NULL, bw.type=c('dist','knn','nen'), D=NULL, resid.type=c('deviance','pearson'), na.action=c(na.omit, na.fail, na.pass), contrasts=NULL) {    
    #Coordinates could be supplied as a character vector or an expression
    coords.is.char = FALSE
    try(coords.is.char <- is.character(coords), silent=TRUE)
    if (coords.is.char) vars = unlist(c("formula", "data", "weights", "na.action"), coords)
    else vars = c("formula", "data", "weights", "na.action")
    
    #Parse the function call, extracting variables, weights, coordinates, and na.action.
    mf = match.call(expand.dots=FALSE)
    m <- match(c("formula", "data", "weights", "na.action"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf[[1]] <- as.name("model.frame")
    mf <- eval.parent(mf)
    mt <- attr(mf, "terms") 
    
    #If the data was provided as a spatial data frame, then extract both the data and the coordinates.
    if (is(data, "Spatial")) {
        if (!missing(coords)) 
            warning("data is Spatial* object, ignoring coords argument")
        coords <- coordinates(data)
        if ((is.null(longlat) || !is.logical(longlat)) && !is.na(is.projected(data)) && !is.projected(data))
            longlat <- TRUE
        else longlat <- FALSE
    } else {
        #Make sure coordinates were specified
        if (missing(coords)) 
            stop("Observation coordinates have to be given")
        
        #Get the coords from the data:
        if (coords.is.char) coords = data[,coords]
        else {
            coords.expression = substitute(coords)
            coords.expression[[1]] = as.name('cbind')
            coords = eval(coords.expression, data)
        }
        
        #Only interpret the coordinates as latitude/longitude values if the longlat variable is TRUE
        if (is.null(longlat) || !is.logical(longlat)) 
            longlat <- FALSE
    }
    
    #Get the data and the weights
    y <- model.response(mf, "numeric")
    x <- model.matrix(mt, mf, contrasts)
    w <- model.weights(mf)
        
    #Get the matrices of distances and weights
    D.coords = rbind(coords, fit.loc)
    if (is.null(D)) {
        n = nrow(D.coords)
        if (longlat) {
            D = as.matrix(earth.dist(D.coords),n,n)
        } else {
            Xmat = matrix(rep(D.coords[,1], times=n), n, n)
            Ymat = matrix(rep(D.coords[,2], times=n), n, n)
            D = sqrt((Xmat-t(Xmat))**2 + (Ymat-t(Ymat))**2)
        }
    }
    
    #Check for problems with the (prior) weights
    if (!is.null(w) && !is.numeric(w)) 
        stop("'weights' must be a numeric vector")
    if (is.null(w)) 
        w <- rep(1, nrow(data))
    if (any(is.na(w))) 
        stop("NAs in weights")
    if (any(w < 0)) 
        stop("negative weights")
    
    #Set some variables that determine how we fit the model
    resid.type = match.arg(resid.type)
    bw.type = match.arg(bw.type)
    varselect.method = match.arg(varselect.method)
    
    #Find the possible range of bandwidths (for use with the adaptive bandwith methods - knn or nen)
    bbox <- cbind(range(coords[, 1]), range(coords[, 2]))
    difmin <- spDistsN1(bbox, bbox[2, ], longlat)[1]
    if (any(!is.finite(difmin))) 
        difmin[which(!is.finite(difmin))] <- 0
    min.dist = difmin / 300
    max.dist = 10 * difmin
    
    #Fit the model:
    res = list()
    res[['model']] = lagr.dispatch(
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
        D=D,
        longlat=longlat,
        varselect.method=varselect.method,
        verbose=verbose,
        bw=bw,
        bw.type=bw.type,
        kernel=kernel,
        min.dist=min.dist,
        max.dist=max.dist,
        tol.loc=tol.loc,
        resid.type=resid.type
    )
    
    #Store results from model fitting:
    if (!tuning) {
        res[['data']] = data
        res[['response']] = as.character(formula[[2]])
        res[['family']] = family
        res[['weights']] = w
        res[['coords']] = coords
        res[['fit.locs']] = fit.loc
        res[['longlat']] = longlat
        res[['kernel']] = kernel
        res[['bw']] = bw
        res[['bw.type']] = bw.type
        res[['varselect.method']] = varselect.method
        
        res[['na.action']] <- attr(mf, "na.action")
        res[['contrasts']] <- attr(x, "contrasts")
        res[['xlevels']] <- .getXlevels(mt, mf)
        #res[['call']] <- cl
        res[['terms']] <- mt
    }
    class(res) = "lagr"
    
    res
}
