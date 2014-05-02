\name{lagr.sel}
\alias{lagr.sel}
\title{Estimate the bandwidth parameter for a LAGR model}
\usage{
lagr.sel(formula, data = list(), family, range = NULL, weights = NULL,
  coords, oracle = NULL, kernel = gwr.Gauss, bw.type = c("dist", "knn",
  "nen"), varselect.method = c("AIC", "BIC", "AICc"), verbose = FALSE,
  longlat = FALSE, tol.loc = .Machine$double.eps^0.25,
  tol.bw = .Machine$double.eps^0.25, parallel = FALSE, interact = FALSE,
  bwselect.method = c("AICc", "GCV", "BICg"), resid.type = c("deviance",
  "pearson"))
}
\description{
Estimate the bandwidth parameter for a LAGR model
}

