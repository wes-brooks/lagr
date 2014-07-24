\name{lagr.tune}
\alias{lagr.tune}
\title{Estimate the bandwidth parameter for a lagr model}
\usage{
lagr.tune(formula, data = list(), family, range = NULL, weights = NULL,
  coords, oracle = NULL, kernel = NULL, bw.type = c("dist", "knn", "nen"),
  varselect.method = c("AIC", "BIC", "AICc"), verbose = FALSE,
  longlat = FALSE, tol.loc = .Machine$double.eps^0.25,
  tol.bw = .Machine$double.eps^0.25, bwselect.method = c("AICc", "GCV",
  "BICg"), resid.type = c("deviance", "pearson"))
}
\arguments{
  \item{formula}{A formula object describing the response
  and the predictor variables.}

  \item{data}{A data frame containing the model-building
  data}

  \item{family}{The exponential family (or Cox model)
  distribution of the response.}

  \item{range}{The allowable range of the bandwidth
  parameter.}

  \item{weights}{Prior weights on the observations. These
  aren't the kernel weights - those will be calculated
  internally.}

  \item{coords}{The coordinates of theobservation
  locations. Rows of \code{coords} must align with rows of
  \code{data}.}

  \item{longlat}{Are the coordinates provided in longitude
  and latitude? Default is \code{FALSE}.}

  \item{kernel}{The kernel function to use for locally
  weighting observations.}

  \item{bw}{The bandwidth to use with the kernel.}

  \item{bw.type}{The type of bandwidth. \code{dist} means
  the bandwidth is a distance, \code{knn} means it is a
  proportion of the \code{n}, and \code{nen} means it is a
  proportion of the global-model sum of squared residuals.}

  \item{varselect.method}{What criterion to minimize during
  variable selection. Options are \code{AIC}, \code{BIC},
  \code{AICc}, and \code{GCV}.}

  \item{bwselect.method}{The name of the bandwidth
  selection criterion (options are \code{AIC}, \code{AICc},
  \code{BIC}, \code{GCV}).}

  \item{resid.type}{What kind of residuals to use. Options
  are \code{pearson} and \code{deviance}.}

  \item{tol.loc}{The allowable error tolerance when
  converting an adaptive bandwidth (\code{knn} or
  \code{nen}) to a distance for each local model.}

  \item{tol.bw}{The allowable error tolerance when finding
  the bandwidth that minimizes the bandwidth selection
  criterion.}
}
\value{
bw The value of bandwidth that minimizes the nbandwidth
selection criterion.

trace A data frame of each bandwidth that was tried during
the optimization, along with the resulting degrees of
freedom used inthe LAGR model and the value of the
bandwidth selection criterion.
}
\description{
\code{lagr.tune} estimates the bandwidth parameter for a
LAGR model.
}
\details{
This method calls \code{lagr} repeatedly via the
\code{optimize} function, searching for the bandwidth that
minimizes a bandwidth selection criterion. It returns the
profiled value of the selection criterion at each bandwidth
that is used in the evaluation.
}

