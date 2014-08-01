#' The Epanechnikov kernel
#' @param R vector of distances to convert to weights
#' @param bw bandwidth of the kernel
#' @return vector of kernel weights, given by \eqn{1 - (R / bw)^2} if \eqn{R < bw}; otherwise 0
#' @export
epanechnikov = function(R, bw) {
    ifelse( R < bw, 1-(R/bw)**2, 0)
}

#' The cubic kernel
#' @param d vector of distances to convert to weights
#' @param bw bandwidth of the kernel
#' @return vector of kernel weights, given by \eqn{1 - (7 * (d / bw)^2 - 8.75 * (d / bw)^3 + 3.5 * (d / bw)^5 - 0.75 * (d / bw)^7} if \eqn{d < bw}; otherwise 0
#' @export
cubic = function(d, bw) {
    ifelse(d<bw, 1 - (7*(d/bw)**2 - 8.75*(d/bw)**3 + 3.5*(d/bw)**5 - 0.75*(d/bw)**7), 0)
}

#' The spherical kernel
#' @param d vector of distances to convert to weights
#' @param bw bandwidth of the kernel
#' @return vector of kernel weights, given by \eqn{1 - 1.5 * (d / bw) + 0.5 * (d / bw)^3} if \eqn{d < bw}; otherwise 0
#' @export
spherical = function(d, bw) {
    ifelse(d<bw, 1 - 1.5*(d/bw) + 0.5*(d/bw)**3, 0)
}

#' The bisquare kernel
#' @param R vector of distances to convert to weights
#' @param bw bandwidth of the kernel
#' @return vector of kernel weights, given by \eqn{(1 - (R / bw)^2)^2} if \eqn{R < bw}; otherwise 0
#' @export
bisquare = function(R, bw) {
    ifelse( R < bw, (1 - (R/bw)**2)**2, 0)
}
