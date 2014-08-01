\name{spherical}
\alias{spherical}
\title{The spherical kernel}
\usage{
spherical(d, bw)
}
\arguments{
  \item{d}{vector of distances to convert to weights}

  \item{bw}{bandwidth of the kernel}
}
\value{
vector of kernel weights, given by \eqn{1 - 1.5 * (d / bw)
+ 0.5 * (d / bw)^3} if \eqn{d < bw}; otherwise 0
}
\description{
The spherical kernel
}

