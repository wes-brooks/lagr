#' Dispatch the fitting of local models to parallel cores, if registered
#' 
#' Loops through estimation locations with foreach, sending each local model to a core for fitting. If zero or one cores are registered, then foreach computes the local models sequentially.
#' 
#' @param x matrix of observed covariates
#' @param y vector of observed responses
#' @param family exponential family distribution of the response
#' @param coords matrix of locations, with each row giving the location at which the corresponding row of data was observed
#' @param fit.loc matrix of locations where the local models should be fitted
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
lagr.dispatch = function(x, y, family, coords, fit.loc, oracle, D, bw, bw.type, verbose, varselect.method, prior.weights, tuning, predict, simulation, kernel, min.bw, max.bw, min.dist, max.dist, tol.loc, lambda.min.ratio, n.lambda, lagr.convergence.tol, lagr.max.iter, resid.type, jacknife=FALSE, bootstrap.index=NULL) {
    if (!is.null(fit.loc)) { coords.fit = fit.loc }
    else { coords.fit = coords }
    n = nrow(coords.fit)

    lagr.object = list()

    #For the adaptive bandwith methods, use a default tolerance if none is specified:
    if (is.null(tol.loc)) {tol.loc = bw / 1000}

    #The knn bandwidth is a proportion of the total prior weight, so compute the total prior weight:
    if (bw.type == 'knn') {
        prior.weights = drop(prior.weights)
        total.weight = sum(prior.weights)
    }

    group.id = attr(x, 'assign')
    
    models = foreach(i=1:n, .errorhandling='pass') %dopar% {
        if (!is.null(fit.loc)) {
            dist = drop(D[nrow(coords)+i,1:nrow(coords)])
        } else { dist = drop(D[i,]) }
        loc = coords.fit[i,]

        #If we are seeking the bandwidth via the jacknife, then remove any observations with zero distance.
        if (jacknife==TRUE || jacknife=='anti') {
            indx = which(dist!=0)
        } else if (jacknife==FALSE) {
            indx = 1:nrow(x)
        } 
        
        if (!is.null(bootstrap.index)) {
            indx = indx[bootstrap.index]
        }
        
        #If this is an anti-jacknife procedure then add the central observation back to the mix:
        if (jacknife=='anti') {
            indx = c(indx, which(dist==0)[1])
        } 
        
        #Use a prespecified distance as the bandwidth?
        if (bw.type == 'dist') {
            bandwidth = bw
            kernel.weights = drop(kernel(dist, bandwidth))

        #Compute the bandwidth that sets the sum of weights around location i equal to bw?
        } else if (bw.type == 'knn') {
            opt = optimize(
                lagr.knn,
                lower=min.dist,
                upper=max.dist,
                maximum=FALSE,
                tol=tol.loc,
                loc=loc,
                coords=coords[indx,],
                kernel=kernel,
                verbose=verbose,
                dist=dist[indx],
                total.weight=total.weight,
                prior.weights=prior.weights[indx],
                target=bw
            )
            bandwidth = opt$minimum
            kernel.weights = drop(kernel(dist, bandwidth))
            
        #Compute the bandwidth so that the sum of the local weighted squared error equals the bw?
        } else if (bw.type == 'nen') {
            opt = optimize(
                lagr.ssr,
                lower=min.dist,
                upper=max.dist, 
                maximum=FALSE,
                tol=tol.loc,
                x=x[indx,],
                y=y[indx],
                group.id=group.id,
                family=family,
                loc=loc,
                coords=coords[indx,],
                dist=dist[indx],
                kernel=kernel,
                target=bw,
                varselect.method=varselect.method,
                resid.type=resid.type,
                oracle=oracle,
                prior.weights=prior.weights[indx],
                verbose=verbose,
                lambda.min.ratio=lambda.min.ratio,
                n.lambda=n.lambda, 
                lagr.convergence.tol=lagr.convergence.tol,
                lagr.max.iter=lagr.max.iter
            )
            bandwidth = opt$minimum
            kernel.weights = drop(kernel(dist, bandwidth))
        }
        
        #If we have specified covariates via an oracle, then use those
        if (!is.null(oracle)) {oracle.loc = oracle[[i]]}
        else {oracle.loc = NULL}

        #Fit the local model
        m = list(tunelist=list('ssr-loc'=list('pearson'=Inf, 'deviance'=Inf), 'df-local'=1), 'sigma2'=0, 'nonzero'=vector(), 'weightsum'=sum(kernel.weights))
        try(m <- lagr.fit.inner(
            x=x[indx,],
            y=y[indx],
            group.id=group.id,
            family=family,
            coords=coords[indx,],
            loc=loc,
            varselect.method=varselect.method,
            tuning=tuning,
            predict=predict,
            simulation=simulation,
            lambda.min.ratio=lambda.min.ratio,
            n.lambda=n.lambda, 
            lagr.convergence.tol=lagr.convergence.tol,
            lagr.max.iter=lagr.max.iter,
            verbose=verbose,
            kernel.weights=kernel.weights[indx],
            prior.weights=prior.weights[indx],
            oracle=oracle.loc)
        )
        if (verbose) {
            cat(paste("For i=", i, "; location=(", paste(round(loc,3), collapse=","), "); bw=", round(bandwidth,3), "; s=", m[['s']], "; dispersion=", round(tail(m[['dispersion']],1),3), "; nonzero=", paste(m[['nonzero']], collapse=","), "; weightsum=", round(m[['weightsum']],3), ".\n", sep=''))
        }
        m[['bw']] = bandwidth
        return(m)
    }
    
    return(models)
}
