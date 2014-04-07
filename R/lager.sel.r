lager.sel = function(formula, data=list(), family, range=NULL, weights=NULL, coords, oracle=NULL, kernel=gwr.Gauss, bw.type=c('dist','knn','nen'), varselect.method=c('AIC','BIC','AICc'), verbose=FALSE, longlat=FALSE, tol.loc=.Machine$double.eps^0.25, tol.bw=.Machine$double.eps^0.25, parallel=FALSE, interact=FALSE, bwselect.method=c('AICc','GCV','BICg'), resid.type=c('deviance','pearson')) {
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
  opt <- optimize(lager.cv.f, interval=c(beta1, beta2), tol=tol.bw, maximum=FALSE,
    formula=formula, coords=coords, env=oo, oracle=oracle, family=family, varselect.method=varselect.method,
    kernel=kernel, verbose=verbose, longlat=longlat, data=data, bw.type=bw.type,
    weights=weights, tol.loc=tol.loc, parallel=parallel, N=1, interact=interact,
    resid.type=resid.type, bwselect.method=bwselect.method)
  trace = oo$trace[!duplicated(oo$trace[,1]),]
  rm(oo)

  bdwt <- opt$minimum
  res <- bdwt
  return(list(bw=res, trace=trace, bwselect.method=bwselect.method, resid.type=resid.type))
}
