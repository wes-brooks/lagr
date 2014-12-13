#' @export
paraboot = function(obj, n.bw, n.draws) {
    if (!is(obj, "lagr.tune")) { stop("Error: obj must be an object of class lagr!") }

    n.loc = length(obj[['fits']])
    response = vector()

    #draw from the EDF interpolated from the bw trace:
    obj$trace = obj$trace[order(obj$trace$bw),]
    ll = -0.5 * obj$trace$loss
    ll = ll - min(ll)
    lik = exp(ll) - exp(min(ll))
    cumLik = cumsum(lik)
    cs = cumsum(c(0,diff(cumLik) * diff(obj$trace$bw) ))

    cumLik = (cs - min(cs)) / diff(range(cs))
    slope = diff(obj$trace$bw) / diff(cumLik)

    len = length(lik)
    obj$trace$bw[1:(len-1)] - slope*lik[1:(len-1)]

#> a[-10,1]-b*a[-10,2] -> int
#> int[5]+b[5]*0.5
#[1] 0.5364039
#> which(diff(0.5<a[,2]))
#Error in which(diff(0.5 < a[, 2])) : argument to 'which' is not logical
#> which(diff(0.5<a[,2])==1)
#[1] 5
#> int[5]+b[5]*0.5


    q = runif(n.bw)
    bws = quantile(obj$trace$b, q, type=4)

    for (i in 1:n.bw) {
        quantile(bw, q, type=4)
    }

    for (i in 1:n.loc) {
        response = c(response, (m[[i]]$coef + (sqrtm(Sigma) %*% y.hat.seed))[1])
    }
}