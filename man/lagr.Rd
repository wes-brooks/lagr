\docType{package}
\name{lagr}
\alias{lagr}
\alias{lagr-package}
\title{Fit a lagr model}
\usage{
lagr(formula, data, family = c("gaussian", "binomial", "poisson", "Cox"),
  weights = NULL, coords, fit.loc = NULL, tuning = FALSE,
  predict = FALSE, simulation = FALSE, oracle = NULL, kernel, bw = NULL,
  varselect.method = c("AIC", "BIC", "AICc"), verbose = FALSE, longlat,
  tol.loc = NULL, bw.type = c("dist", "knn", "nen"), D = NULL,
  resid.type = c("deviance", "pearson"))
}
\arguments{
  \item{formula}{A formula object describing the response
  and the predictor variables.}

  \item{data}{A data frame containing the model-building
  data}

  \item{family}{The exponential family (or Cox model)
  distribution of the response.}

  \item{weights}{Prior weights on the observations. These
  aren't the kernel weights - those will be calculated
  internally.}

  \item{coords}{The coordinates of theobservation
  locations. Rows of \code{coords} must align with rows of
  \code{data}.}

  \item{fit.loc}{(optional) You can specify the locations
  to fit the local models via \code{fit.loc}. Otherwise,
  the locations in \code{coords} will be used.}

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

  \item{resid.type}{What kind of residuals to use. Options
  are \code{pearson} and \code{deviance}.}
}
\value{
models A list containing the local models.
}
\description{
\code{lagr} fits a model via local, adaptive, grouped
regularization.
}
\details{
This method fits a local model at each location indicated
by \code{fit.loc}.
}
\references{
Garrett Grolemund, Hadley Wickham (2011). Dates and Times
Made Easy with lubridate. Journal of Statistical Software,
40(3), 1-25. \url{http://www.jstatsoft.org/v40/i03/}.
}

