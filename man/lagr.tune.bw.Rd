\name{lagr.tune.bw}
\alias{lagr.tune.bw}
\title{Evaluate the bandwidth selection criterion for a given bandwidth}
\usage{
lagr.tune.bw(formula, data, weights, family, bw, kernel, coords, longlat, env,
  oracle, varselect.method, tol.loc, bw.type, bwselect.method, resid.type,
  verbose)
}
\arguments{
  \item{formula}{symbolic representation of the model}

  \item{data}{data frame containing observations of all the
  terms represented in the formula}

  \item{weights}{vector of prior observation weights (due
  to, e.g., overdispersion). Not related to the kernel
  weights.}

  \item{family}{exponential family distribution of the
  response}

  \item{bw}{bandwidth for the kernel}

  \item{kernel}{kernel function for generating the local
  observation weights}

  \item{coords}{matrix of locations, with each row giving
  the location at which the corresponding row of data was
  observed}

  \item{longlat}{\code{TRUE} indicates that the coordinates
  are specified in longitude/latitude, \code{FALSE}
  indicates Cartesian coordinates. Default is
  \code{FALSE}.}

  \item{varselect.method}{criterion to minimize in the
  regularization step of fitting local models - options are
  \code{AIC}, \code{AICc}, \code{BIC}, \code{GCV}}

  \item{tol.loc}{tolerance for the tuning of an adaptive
  bandwidth (e.g. \code{knn} or \code{nen})}

  \item{bw.type}{type of bandwidth - options are
  \code{dist} for distance (the default), \code{knn} for
  nearest neighbors (bandwidth a proportion of \code{n}),
  and \code{nen} for nearest effective neighbors (bandwidth
  a proportion of the sum of squared residuals from a
  global model)}

  \item{bwselect.method}{criterion to minimize when tuning
  bandwidth - options are \code{AICc}, \code{BICg}, and
  \code{GCV}}

  \item{resid.type}{type of residual to use (relevant for
  non-gaussian response) - options are \code{deviance} and
  \code{pearson}}

  \item{verbose}{print detailed information about our
  progress?}
}
\value{
value of the \code{bwselect.method} criterion for the given
bandwidth
}
\description{
Evaluate the bandwidth selection criterion for a given
bandwidth
}

