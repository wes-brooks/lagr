#' Estimate the bandwidth parameter for a lagr model
#'
#' \code{lagr.sel} estimates the bandwidth parameter for a LAGR model.
#'
#' This method calls \code{lagr} repeatedly via the \code{optimize} function, searching for the bandwidth that minimizes a bandwidth selection criterion. It returns the profiled value of the selection criterion at each bandwidth that is used in the evaluation.
#'
#' @param formula A formula object describing the response and the predictor variables.
#' @param data A data frame containing the model-building data
#' @param family The exponential family (or Cox model) distribution of the response.
#' @param range The allowable range of the bandwidth parameter.
#' @param weights Prior weights on the observations. These aren't the kernel weights - those will be calculated internally.
#' @param coords The coordinates of theobservation locations. Rows of \code{coords} must align with rows of \code{data}.
#' @param longlat Are the coordinates provided in longitude and latitude? Default is \code{FALSE}.
#' @param kernel The kernel function to use for locally weighting observations.
#' @param bw The bandwidth to use with the kernel.
#' @param bw.type The type of bandwidth. \code{dist} means the bandwidth is a distance, \code{knn} means it is a proportion of the \code{n}, and \code{nen} means it is a proportion of the global-model sum of squared residuals.
#' @param varselect.method What criterion to minimize during variable selection. Options are \code{AIC}, \code{BIC}, \code{AICc}, and \code{GCV}.
#' @param bwselect.method The name of the bandwidth selection criterion (options are \code{AIC}, \code{AICc}, \code{BIC}, \code{GCV}).
#' @param resid.type What kind of residuals to use. Options are \code{pearson} and \code{deviance}.
#' @param tol.loc The allowable error tolerance when converting an adaptive bandwidth (\code{knn} or \code{nen}) to a distance for each local model.
#' @param tol.bw The allowable error tolerance when finding the bandwidth that minimizes the bandwidth selection criterion.
#' 
#' @return bw The value of bandwidth that minimizes the nbandwidth selection criterion.
#' @return trace A data frame of each bandwidth that was tried during the optimization, along with the resulting degrees of freedom used inthe LAGR model and the value of the bandwidth selection criterion.
#' 
#' @export
lagr.sel = function(formula, data=list(), family, range=NULL, weights=NULL, coords, oracle=NULL, kernel=NULL, bw.type=c('dist','knn','nen'), varselect.method=c('AIC','BIC','AICc'), verbose=FALSE, longlat=FALSE, tol.loc=.Machine$double.eps^0.25, tol.bw=.Machine$double.eps^0.25, bwselect.method=c('AICc','GCV','BICg'), resid.type=c('deviance','pearson')) {
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
  opt <- optimize(lagr.cv.f, interval=c(beta1, beta2), tol=tol.bw, maximum=FALSE,
    formula=formula, coords=coords, env=oo, oracle=oracle, family=family, varselect.method=varselect.method,
    kernel=kernel, verbose=verbose, longlat=longlat, data=data, bw.type=bw.type,
    weights=weights, tol.loc=tol.loc,
    resid.type=resid.type, bwselect.method=bwselect.method)
  trace = oo$trace[!duplicated(oo$trace[,1]),]
  rm(oo)

  bdwt <- opt$minimum
  res <- bdwt
  return(list(bw=res, trace=trace, bwselect.method=bwselect.method, resid.type=resid.type))
}
