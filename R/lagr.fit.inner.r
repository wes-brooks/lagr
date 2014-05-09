#' Fit a LAGR model
lagr.fit.inner = function(x, y, coords, loc, family, varselect.method, tuning, predict, simulation, verbose, kernel.weights=NULL, prior.weights=NULL, longlat=FALSE, N=1, oracle=oracle) {

    #Find which observations were made at the model location  
    colocated = which(round(coords[,1],5) == round(as.numeric(loc[1]),5) & round(coords[,2],5) == round(as.numeric(loc[2]),5))

    #Bail now if only one observation has weight:
    if (sum(kernel.weights)==length(colocated)) {
        return(list('tunelist'=list('df-local'=1, 'ssr-loc'=list('pearson'=Inf, 'deviance'=Inf))))
    }

    #Establish groups for the group lasso
    vargroup = 1:ncol(x)
  
    #Compute the covariate-by-location interactions
    raw.names = colnames(x)
    interact.names = vector()
    for (l in 1:length(raw.names)) {
        interact.names = c(interact.names, paste(raw.names[l], ":", colnames(coords)[1], sep=""))
        interact.names = c(interact.names, paste(raw.names[l], ":", colnames(coords)[2], sep=""))
    }

    interacted = matrix(0, ncol=2*ncol(x), nrow=nrow(x))
    for (k in 1:ncol(x)) {
        interacted[,2*(k-1)+1] = x[,k]*(coords[,1]-loc[[1]])
        interacted[,2*k] = x[,k]*(coords[,2]-loc[[2]])
        vargroup = c(vargroup, vargroup[k], vargroup[k])
    }

    x.interacted = cbind(x, interacted)
    colnames(x.interacted) = c(raw.names, interact.names)

    #Combine prior weights and kernel weights
    w <- prior.weights * kernel.weights
    weighted = which(w>0)
    n.weighted = length(weighted)
  
    #Limit our attention to the observations with nonzero weight
    xx = as.matrix(x[weighted,])
    xx.interacted = as.matrix(x.interacted[weighted,])
    yy = as.matrix(y[weighted])
    w = w[weighted]

    #Instantiate objects to store our output
    int.list = list()
    coef.list = list()
    tunelist = list()

    #Loop through N bootstrap replicates
    for (i in 1:N) {
        if (i==N) {
            #Final permutation is the original ordering of the data:
            permutation = 1:n.weighted
        } else {
            #Generate permutations to preserve the total local weight:
            tot.w = sum(w)
            permutation = vector()
            while (sum(w[permutation]) < tot.w) {
                permutation = c(permutation, sample(1:n.weighted, size=1))
            }

            permutation = permutation[1:which.min(tot.w - cumsum(w[permutation]))]
        }

        #Sample the data according to the permutation:
        colocated = which(kernel.weights[weighted][permutation]==1)
        xxx = xx.interacted[permutation,]
        yyy = yy[permutation]
        sumw = sum(w[permutation])
        
        #Use the adaptive group lasso to produce a local model:
        model = SGL(data=list(x=xxx, y=yyy), weights=w[permutation], index=vargroup, maxit=100, standardize=FALSE, alpha=0, delta=2, nlam=20, min.frac=0.00001, thresh=0.01, adaptive=TRUE)
    
        vars = apply(as.matrix(model[['beta']]), 2, function(x) {which(x!=0)})
        df = model[['results']][['df']] + 1 #Add one because we must estimate the scale parameter.

        if (sumw > ncol(x)) {
            #Extract the fitted values for each lambda:
            fitted = model[['results']][['fitted']]
            s2 = sum(w[permutation]*(model[['results']][['residuals']][,ncol(fitted)])**2) / (sumw - df) 
      
            #Compute the loss (varies by family)
            #loss = model[[varselect.method]]
            if (varselect.method == 'AIC') {penalty = 2*df}
            if (varselect.method == 'AICc') {penalty = 2*df + 2*df*(df+1)/(sumw - df - 1)}
            if (varselect.method == 'BIC') {penalty = sumw*df}


            #Assuming scale from the largest model:
            #loss = sumw * log(s2) + apply(model[['results']][['residuals']], 2, function(x) sum(w[permutation]*x**2))/s2 + penalty

            #Estimating scale in penalty formula:
            loss = sumw * (log(apply(model[['results']][['residuals']], 2, function(x) sum(w[permutation]*x**2))) - log(sumw) + 1) + penalty

            #Estimating the loss only at the modeling location (not the total local loss)
            #loss = (model[['results']][['residuals']][colocated,])**2 + penalty*w[permutation][colocated] / sumw

            #Pick the lambda that minimizes the loss:
            k = which.min(loss)
            fitted = fitted[,k]
            localfit = fitted[colocated]
            df = df[k]
            if (k > 1) {
                varset = vars[[k]]
            } else {
                varset = NULL
            }
      
            if (length(colocated)>0) {
                tunelist[['ssr-loc']] = list()
                tunelist[['ssr']] = list()
        
                #Pearson residuals:
                if (family=='gaussian') {
                    tunelist[['ssr-loc']][['pearson']] = sum((w[permutation]*(fitted - yyy)**2)[colocated])
                    tunelist[['ssr']][['pearson']] = sum(w[permutation]*(fitted - yyy)**2)
                } else if (family=='poisson') {
                    tunelist[['ssr-loc']][['pearson']] = sum((w[permutation]*(yyy - fitted)**2/fitted)[colocated])
                    tunelist[['ssr']][['pearson']] = sum(w[permutation]*(fitted - yyy)**2/fitted)
                } else if (family=='binomial') {
                    tunelist[['ssr-loc']][['pearson']] = sum((w[permutation]*(yyy - fitted)**2/(fitted*(1-fitted)))[colocated])
                    tunelist[['ssr']][['pearson']] = sum(w[permutation]*(fitted - yyy)**2/(fitted*(1-fitted)))
                }

                #Deviance residuals:
                if (family=='gaussian') {
                    tunelist[['ssr-loc']][['deviance']] = sum((w[permutation]*(fitted - yyy)**2)[colocated])
                    tunelist[['ssr']][['deviance']] = sum(w[permutation]*(fitted - yyy)**2)
                } else if (family=='poisson') {
                    tunelist[['ssr-loc']][['deviance']] = sum((2*w[permutation]*(ylogy(yyy) - yyy*log(fitted) - (yyy-fitted)))[colocated])
                    tunelist[['ssr']][['deviance']] = sum(2*w[permutation]*(ylogy(yyy) - yyy*log(fitted) - (yyy-fitted)))
                } else if (family=='binomial') {
                    tunelist[['ssr-loc']][['deviance']] = sum((2*w[permutation]*(ylogy(yyy) - yyy*log(fitted) - ylogy(1-yyy) + (1-yyy)*log(1-fitted)))[colocated])
                    tunelist[['ssr']][['deviance']] = sum(2*w[permutation]*(ylogy(yyy) - yyy*log(fitted) - ylogy(1-yyy) + (1-yyy)*log(1-fitted)))
                }

                #Compute the dispersion parameter:
                if (family=='gaussian') { tunelist[['s2']] = s2 }
                else if (family=='poisson') { tunelist[['s2']] = summary(m)$dispersion }
                else if (family=='binomial') { tunelist[['s2']] = 1 }
        
                #Prepare some outputs for the bandwidth-finding scheme:
                tunelist[['n']] = sumw
                tunelist[['df']] = df
                tunelist[['df-local']] = df*w[permutation][colocated] / sumw
            } else {
                loss.local = NA
            }                   
        } else {
            fitted = rep(meany, length(permutation))
            s2 = 0
            loss = Inf
            loss.local = c(Inf)   
            localfit = meany
        }
    
        #Get the coefficients:
        coefs = t(rbind(model[['intercept']], model[['beta']]))[k,]
        coefs = Matrix(coefs, ncol=1)
        rownames(coefs) = c("(Intercept)", colnames(xxx))
        coefs = coefs[c("(Intercept)", raw.names),]

        coef.list[[i]] = coefs
    }
  
    if (tuning) {
        return(list(tunelist=tunelist, s=k, sigma2=s2, nonzero=raw.names[unique(vargroup[vars[[k]]])], weightsum=sumw, loss=loss))
    } else if (predict) {
        return(list(tunelist=tunelist, coef=coefs, weightsum=sumw, s=k, sigma2=s2, nonzero=raw.names[unique(vargroup[vars[[k]]])]))
    } else if (simulation) {
        return(list(tunelist=tunelist, coef=coefs, coeflist=coef.list, s=k, sigma2=s2, fitted=localfit, nonzero=raw.names[unique(vargroup[vars[[k]]])], actual=yyy[colocated], weightsum=sumw, loss=loss))
    } else {
        return(list(model=model, loss=loss, coef=coefs, coeflist=coef.list, nonzero=raw.names[unique(vargroup[vars[[k]]])], s=k, loc=loc, df=df, loss.local=loss, sigma2=s2, N=N, fitted=localfit, weightsum=sumw))
    }
}
