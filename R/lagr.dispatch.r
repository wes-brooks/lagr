lagr.dispatch = function(x, y, family, coords, fit.loc, longlat, oracle, D, bw.type, verbose, varselect.method, prior.weights, tuning, predict, simulation, kernel, target, min.bw, max.bw, min.dist, max.dist, tol.loc, N, interact, resid.type) {
  if (!is.null(fit.loc)) { coords.unique = unique(fit.loc) }
  else { coords.unique = unique(coords) }
  n = dim(coords.unique)[1]
  
  lagr.object = list()
  
  #For the adaptive bandwith methods, use a default tolerance if none is specified:
  if (is.null(tol.loc)) {tol.loc = target / 1000}
  
  #The knn bandwidth is a proportion of the total prior weight, so compute the total prior weight:
  if (bw.type == 'knn') {
    prior.weights = drop(prior.weights)
    max.weights = rep(1, length(prior.weights))
    total.weight = sum(max.weights * prior.weights)
  }
  
  models = foreach(i=1:n, .packages=c('SGL'), .errorhandling='remove') %dopar% {
    loc = coords.unique[i,]
    dist = drop(D[i,])
    
    #Use a prespecified distance as the bandwidth?
    if (bw.type == 'dist') {
      bandwidth = target
      kernel.weights = drop(kernel(dist, bandwidth))
      
    #Compute the bandwidth that sets the sum of weights around location i equal to target?
    } else if (bw.type == 'knn') {
      opt = optimize(lagr.knn,
                     lower=min.dist,
                     upper=max.dist,
                     maximum=FALSE,
                     tol=tol.loc,
                     coords=coords,
                     loc=loc,
                     kernel=kernel,
                     verbose=verbose,
                     dist=dist,
                     total.weight=total.weight,
                     prior.weights=prior.weights,
                     target=target)
      bandwidth = opt$minimum
      kernel.weights = drop(kernel(dist, bandwidth))
      
    #Compute the bandwidth so that the sum of the local weighted squared error equals the target?
    } else if (bw.type == 'nen') {
      opt = optimize(lagr.ssr,
                     lower=min.dist,
                     upper=max.dist, 
                     maximum=FALSE,
                     tol=tol.loc,
                     x=x,
                     y=y,
                     coords=coords,
                     loc=loc,
                     resid.type=resid.type,
                     kernel=kernel,
                     verbose=verbose,
                     dist=dist,
                     varselect.method=varselect.method,
                     family=family,
                     oracle=oracle,
                     prior.weights=prior.weights,
                     target=target,
                     interact=interact)
      bandwidth = opt$minimum
      kernel.weights = drop(kernel(dist, bandwdth))
    }
    
    #Fit the local model
    m = lagr.fit.inner(x=x,
                        y=y,
                        family=family,
                        coords=coords,
                        loc=loc,
                        N=N,
                        varselect.method=varselect.method,
                        tuning=tuning,
                        predict=predict,
                        simulation=simulation,
                        verbose=verbose,
                        kernel.weights=kernel.weights,
                        prior.weights=prior.weights,
                        interact=interact, 
                        oracle=oracle)
    
    if (verbose) {
      cat(paste("For i=", i, "; location=(", paste(round(loc,3), collapse=","), "); bw=", round(bandwidth,3), "; s=", m[['s']], "; sigma2=", round(tail(m[['sigma2']],1),3), "; nonzero=", paste(m[['nonzero']], collapse=","), "; weightsum=", round(m[['weightsum']],3), ".\n", sep=''))
    }
    
    return(m)
  }
  
  lagr.object[['models']] = models
  
  if (tuning) { }
  else if (predict) { }
  else if (simulation) { }
  else {lagr.object[['coords']] = coords}
  
  class(lagr.object) = 'lagr.object'
  return(lagr.object)
}
