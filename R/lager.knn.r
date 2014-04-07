lager.knn = function(bw, coords, indx=NULL, loc, dist, verbose, prior.weights, total.weight, kernel, target) {
    kernel.weights = kernel(dist, bw)

    if (!is.null(indx)) {
        kernel.weights = kernel.weights[indx]
    }

    w = kernel.weights * prior.weights
    prop = sum(w)/total.weight

    return(abs(prop-target))
}
