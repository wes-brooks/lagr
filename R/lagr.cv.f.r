lagr.cv.f = function(formula, data, weights, family, bw, coords, kernel, env, oracle, varselect.method, verbose, longlat, tol.loc, bw.type, bwselect.method, resid.type) {    
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
                      resid.type=resid.type)
  
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
