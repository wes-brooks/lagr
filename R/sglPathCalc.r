betterPathCalc <- function(data, index, weights, adaweights, min.frac=0.05, nlam=20, type="linear") {
    reset <- 10
    step <- 1
    gamma <- 0.8
    
    inner.iter <- 1000
    outer.iter <- 1000
    thresh = 10^(-3)
    outer.thresh = thresh
    
    n <- nrow(data$x)
    if (type=="linear") {
        X <- data$x
        resp <- data$y
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
        num.groups <- length(groups)
        range.group.ind <- rep(0,(num.groups+1))
        for (i in 1:num.groups) {
            range.group.ind[i] <- min(which(index == groups[i])) - 1
        }
        range.group.ind[num.groups + 1] <- ncol(X)
        group.length <- diff(range.group.ind)
    }
    
    if (type=="logit") {
        X <- data$x
        y <- data$y
        n <- nrow(X)
        p <- ncol(X)
        
        ## Setting up group lasso stuff ## 
        ord <- order(index)
        index <- index[ord]
        X <- X[,ord]
        unOrd <- match(1:length(ord),ord)
        
        ## Coming up with other C++ info ##
        groups <- unique(index)
        num.groups <- length(groups)
        range.group.ind <- rep(0,(num.groups+1))
        for (i in 1:num.groups) {
            range.group.ind[i] <- min(which(index == groups[i])) - 1
        }
        range.group.ind[num.groups + 1] <- ncol(X)
        
        group.length <- diff(range.group.ind)
        beta.naught <- rep(0,ncol(X))
        beta <- beta.naught
        
        beta.is.zero <- rep(1, num.groups)
        beta.old <- rep(0, ncol(X))
        betas <- matrix(0, nrow = ncol(X), ncol = nlam)
        
        eta <- rep(0,n)
        intercepts <- mean(y)
        eta = eta + intercepts
        m.y <- mean(y)
        resp <- m.y*m.y*(1-m.y) - (y-m.y)
    }
    
    if (type=="cox") {
        covariates <- data$x
        n <- nrow(covariates)
        p <- ncol(covariates)  
        time <- data$time
        status <- data$status
        
        ## Ordering Response and Removing any Censored obs before first death ##
        death.order <- order(time)
        ordered.time <- sort(time)  
        
        X <- covariates[death.order,]  
        ordered.status <- status[death.order]
        
        first.blood <- min(which(ordered.status == 1))
        
        X <- X[first.blood:n,]
        ordered.status <- ordered.status[first.blood:n]
        ordered.time <- ordered.time[first.blood:n]
        death.order <- death.order[first.blood:n]
        n <- n-first.blood+1
        
        death.times <- unique(ordered.time[which(ordered.status == 1)])  ## Increasing list of times when someone died (censored ends not included) ##
        
        ## Calculating Risk Sets ##
        risk.set <- rep(0,n)
        for (i in 1:n) {
            risk.set[i] <- max(which(death.times <= ordered.time[i]))
        }
        
        ## Calculating risk set beginning/ending indices ##
        risk.set.ind <- rep(0,(length(death.times)+1))  
        for(i in 1:length(death.times)){
            risk.set.ind[i] <- min(which(ordered.time >= death.times[i]))
        }
        risk.set.ind[length(risk.set.ind)] <- length(ordered.time) + 1
        
        ## Calculating number of deaths at each death time ##
        num.deaths <- rep(0,length(death.times))
        for (i in 1:length(ordered.time)) {
            if (ordered.status[i] == 1) {
                num.deaths[which(death.times == ordered.time[i])] <-  num.deaths[which(death.times == ordered.time[i])] + 1
            }
        }
        
        ## Finding death indices and number of deaths ##
        death.index <- which(ordered.status == 1)
        total.deaths <- length(death.index)
        
        ## Setting up group lasso stuff ##  
        ord <- order(index)
        index <- index[ord]
        X <- X[,ord]
        unOrd <- match(1:length(ord),ord)
        
        ## Coming up with other C++ info ## 
        groups <- unique(index)
        num.groups <- length(groups)
        range.group.ind <- rep(0,(num.groups+1))
        for (i in 1:num.groups) {
            range.group.ind[i] <- min(which(index == groups[i])) - 1
        }
        range.group.ind[num.groups + 1] <- ncol(X)
        
        group.length <- diff(range.group.ind)
        beta.naught <- rep(0,ncol(X))
        beta <- beta.naught
        
        beta.is.zero <- rep(1, num.groups)
        beta.old <- rep(0, ncol(X))
        beta <- array(0, c(ncol(X),nlam,nlam))
        
        beta.is.zero <- rep(1, num.groups)
        
        eta <- rep(0,n)
        
        ## DONE SETTING UP COX MODEL STUFF
        junk1 <- .C("Cox", riskSetInd = as.integer(risk.set.ind), riskSet = as.integer(risk.set), numDeath = as.integer(num.deaths), status = as.integer(ordered.status), ndeath = as.integer(length(death.times)), nrow = as.integer(n), ncol = as.integer(p), beta = as.double(rep(0,p)), eta = as.double(rep(0,n)), y = as.double(rep(0,n)), weights = as.double(rep(0,n)))
        
        resp <- junk1$y * junk1$weights
    }
    
    lambda.max <- rep(0,num.groups)
    
    #Pure group lasso:
    for (i in 1:num.groups) {
        if (adaweights[i] > 0) {
            ind <- groups[i]
            X.fit <- X[,which(index == ind)]
            cors <- t(X.fit) %*% diag(weights) %*% resp / adaweights[i]
            lambda.max[i] <- sqrt(sum(cors^2)) / sqrt(group.length[i])
        }
    }
    
    max.lam <- max(lambda.max) * 1.1 #The factor of 1.1 is to make sure the first lambda has only unpenalized covariates.
    min.lam <- min.frac * max.lam
    lambdas <- exp(seq(log(max.lam), log(min.lam), length.out=nlam))
    return(lambdas/sum(weights))
}
