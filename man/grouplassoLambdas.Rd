\name{grouplassoLambdas}
\alias{grouplassoLambdas}
\title{Compute the lambdas (adaptive group lasso tuning parameters) at which the local model will be fitted.}
\usage{
grouplassoLambdas(data, index, weights, adaweights, min.frac = 0.05,
  nlam = 20, type = "linear")
}
\arguments{
  \item{data}{list consisting of \code{x}, a matrix of
  grouped lcoal covariates, and \code{y}, the vector of
  observed responses}

  \item{index}{vector indicating group membership for each
  column of \code{x}}

  \item{weights}{vector of observations weights}

  \item{adaweights}{vector of adaptive group penalty
  multipliers (computed from the unpenalized coefficients)}

  \item{min.frac}{ratio of \code{lambda.max} to
  \code{lambda.min}}

  \item{nlam}{how many lambdas to compute along the
  sequence from \code{lambda.max} to \code{lambda.min}}

  \item{type}{string indicating the kind of regression
  model, defaults to \code{"linear"}}
}
\value{
vector of adaptive group lasso tuning parameters at which
to compute the local LAGR coefficients
}
\description{
Finds \code{lambda.max}, the smallest lambda that shrinks
all the penalized groups to zero. The rest of the lambdas
form a log-linear sequence from \code{lambda.max} to
\code{lambda.min = min.frac*lambda.max}.
}

