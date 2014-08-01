\name{epanechnikov}
\alias{epanechnikov}
\title{The Epanechnikov kernel}
\usage{
epanechnikov(R, bw)
}
\arguments{
  \item{R}{vector of distances to convert to weights}

  \item{bw}{bandwidth of the kernel}
}
\value{
vector of kernel weights, given by \eqn{1 - (R / bw)^2} if
\eqn{R < bw}; otherwise 0
}
\description{
The Epanechnikov kernel
}

