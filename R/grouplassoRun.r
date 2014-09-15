#' Get the sequence of tuning parameters and then make the call to the C code that implements the adaptive group lasso
#' 
grouplassoRun <-
    function(data, index, weights, adaweights, family, thresh=0.0001, nlam=20, lambda=NULL, inner.iter=100, outer.iter=100, outer.thresh=0.0001, optim.step.size=0.8, reset=10, min.frac=0.05, verbose=FALSE) {
        lambdaSeek <- grouplassoLambdas(data=data, index=index, weights=weights, family=family, min.frac=min.frac, nlam=nlam, type="linear", adaweights=adaweights)
        if (is.null(lambda)) {
            lambda = lambdaSeek$lambda
        }
        
        X <- data$x
        y <- as.vector(data$y)
        n <- nrow(X)
        p <- ncol(X)
        
        #put the groups and adaweights in numerical order
        groups <- unique(index)
        ord.g = order(groups)
        groups = groups[ord.g]
        adaweights = adaweights[ord.g]
        
        #Reorder columns of X so that groups are contiguous
        ord <- order(index)
        index <- index[ord]
        X <- X[,ord]
        unOrd <- match(1:length(ord),ord)
        
        ## Coming up with other C++ info ##
        inner.iter = as.integer(inner.iter)
        outer.iter = as.integer(outer.iter)
        n = as.integer(nrow(X))
        p = as.integer(ncol(X))
        reset = as.integer(reset)
        
        num.groups = as.integer(length(groups))
        range.group.ind = rep(0,(num.groups+1))
        for(i in 1:num.groups){
            range.group.ind[i] <- min(which(index == groups[i])) - 1
        }
        range.group.ind[num.groups+1] <- ncol(X)
        range.group.ind = as.integer(range.group.ind)
        
        group.length <- as.integer(diff(range.group.ind))
        
        ## DONE SETTING UP C STUFF ##
        
        nlam = length(lambda)
        beta <- lambdaSeek$beta.unpen
        eta <- lambdaSeek$eta.unpen
        
        junk <- rcppLinNest(X = X,
                            y = y,
                            w = weights,
                            linkinv=family$linkinv,
                            devfun=family$dev.resids,
                            adaweights = adaweights,
                            nrow = n,
                            ncol = p,
                            numGroup = num.groups,
                            rangeGroupInd = range.group.ind,
                            groupLen = group.length,
                            lambda = lambda,
                            beta = beta,
                            innerIter = inner.iter,
                            outerIter = outer.iter,
                            thresh = thresh,
                            outerThresh = outer.thresh,
                            eta = eta,
                            gamma = optim.step.size,
                            betaIsZero = as.integer(rep(1, num.groups)),
                            reset = reset
        )
        
        return(list(beta=beta[unOrd,], lambda=lambda))
    }
