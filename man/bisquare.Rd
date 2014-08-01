\name{bisquare}
\alias{bisquare}
\title{The bisquare kernel}
\usage{
bisquare(R, bw)
}
\arguments{
  \item{R}{vector of distances to convert to weights}

  \item{bw}{bandwidth of the kernel}
}
\value{
vector of kernel weights, given by \eqn{(1 - (R / bw)^2)^2}
if \eqn{R < bw}; otherwise 0
}
\description{
The bisquare kernel
}

