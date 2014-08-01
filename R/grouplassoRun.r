#' Get the sequence of tuning parameters and then make the call to the C code that implements the adaptive group lasso
#' 
grouplassoRun <-
    function(data, index, weights, adaweights, thresh=0.0001, nlam=20, lambda=NULL, inner.iter=100, outer.iter=100, outer.thresh=0.0001, gamma=0.8, momentum=1, reset=10, min.frac=0.05, verbose=FALSE) {
        if (is.null(lambda)) {
            lambda <- grouplassoLambdas(data=data, index=index, weights=weights, min.frac=min.frac, nlam=nlam, type="linear", adaweights=adaweights)
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
        beta.old <- rep(0, ncol(X))
        beta.is.zero <- rep(1, num.groups)
        beta <- matrix(0, nrow=ncol(X), ncol=nlam)
        
        beta.is.zero <- as.integer(rep(1, num.groups))
        beta.old <- rep(0, ncol(X))
        eta <- rep(0,n)
        
        junk <- rcppLinNest(X = X,
                            y = y,
                            w = weights,
                            link=identityLinkCpp,
                            loglik=linLogLik,
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
                            gamma = gamma,
                            betaIsZero = beta.is.zero,
                            momentum = momentum,
                            reset = reset
        )
        
        return(list(beta=beta[unOrd,], lambda=lambda))
    }
