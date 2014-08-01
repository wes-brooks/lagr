\name{lagr.dispatch}
\alias{lagr.dispatch}
\title{Dispatch the fitting of local models to parallel cores, if registered}
\usage{
lagr.dispatch(x, y, family, coords, fit.loc, longlat, oracle, D, bw.type,
  verbose, varselect.method, prior.weights, tuning, predict, simulation, kernel,
  target, min.bw, max.bw, min.dist, max.dist, tol.loc, resid.type)
}
\arguments{
  \item{x}{matrix of observed covariates}

  \item{y}{vector of observed responses}

  \item{family}{exponential family distribution of the
  response}

  \item{coords}{matrix of locations, with each row giving
  the location at which the corresponding row of data was
  observed}

  \item{fit.loc}{matrix of locations where the local models
  should be fitted}

  \item{longlat}{\code{TRUE} indicates that the coordinates
  are specified in longitude/latitude, \code{FALSE}
  indicates Cartesian coordinates. Default is
  \code{FALSE}.}

  \item{kernel}{kernel function for generating the local
  observation weights}

  \item{bw}{bandwidth parameter}

  \item{bw.type}{type of bandwidth - options are
  \code{dist} for distance (the default), \code{knn} for
  nearest neighbors (bandwidth a proportion of \code{n}),
  and \code{nen} for nearest effective neighbors (bandwidth
  a proportion of the sum of squared residuals from a
  global model)}

  \item{tol.loc}{tolerance for the tuning of an adaptive
  bandwidth (e.g. \code{knn} or \code{nen})}

  \item{varselect.method}{criterion to minimize in the
  regularization step of fitting local models - options are
  \code{AIC}, \code{AICc}, \code{BIC}, \code{GCV}}

  \item{resid.type}{type of residual to use (relevant for
  non-gaussian response) - options are \code{deviance} and
  \code{pearson}}

  \item{tuning}{logical indicating whether this model will
  be used to tune the bandwidth, in which case only the
  tuning criteria are returned}

  \item{D}{pre-specified matrix of distances between
  locations}

  \item{verbose}{print detailed information about our
  progress?}
}
\description{
Loops through estimation locations with foreach, sending
each local model to a core for fitting. If zero or one
cores are registered, then foreach computes the local
models sequentially.
}

