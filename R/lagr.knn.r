#' Calculate the proportion of observations that lie within a given bandwidth
#' 
#' Facilitates a \code{knn}-type bandwidth by calculating the proportion \eqn{sum(w) / n}. Called repeatedly by \code{optimize} in \code{lagr.dispatch} to home in on the desired \code{knn}.
#' 
#' @param bw bandwidth for which to compute \code{sum(w)}
#' @param loc location around which to center the kernel
#' @param coords matrix of observation locations
#' @param dist vector of distances from central location to the observation locations
#' @param kernel kernel function for generating the local observation weights
#' @param target targeted \code{knn} bandwidth
#' @param prior.weights vector of prior observation weights provided by the user
#' @param total.weight sum of prior weights
#' @param verbose print detailed information about our progress?
#' 
#' @return difference between the calculated \code{sum(w) / total.weight} and the target
#' 
lagr.knn = function(bw, loc, coords, dist, kernel, target, prior.weights, total.weight, verbose) {
    kernel.weights = kernel(dist, bw)
    w = kernel.weights * prior.weights
    prop = sum(w)/total.weight

    #return the difference between the calculated sum(w) and the target
    return(abs(prop-target))
}
