lagr.ssr = function(bw, x, y, family, coords, loc, dist, verbose, prior.weights, kernel, target, varselect.method, interact, oracle, resid.type) {
  #Calculate the local weights:
  kernel.weights = drop(kernel(dist, bw))
  
  lagr.object = lagr.fit.inner(x=x,
                                 y=y,
                                 family=family,
                                 coords=coords,
                                 loc=loc,
                                 varselect.method=varselect.method,
                                 interact=interact,
                                 predict=TRUE,
                                 tuning=FALSE,
                                 simulation=FALSE,
                                 verbose=verbose,
                                 kernel.weights=kernel.weights,
                                 prior.weights=prior.weights,
                                 oracle=oracle)
  
  loss = lagr.object[['ssr']][[resid.type]]
  if (verbose) { cat(paste('loc:(', paste(round(loc,3), collapse=","), '), target: ', round(target,3), ', bw:', round(bw,3), ', ssr:', round(loss,3), ', miss:', round(abs(loss-target),3), '\n', sep="")) }

  return(abs(loss-target))
}
