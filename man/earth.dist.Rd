\name{earth.dist}
\alias{earth.dist}
\title{Find the distance between locations specified by latitude/longitude in degrees.}
\usage{
earth.dist(locs, dist = TRUE)
}
\arguments{
  \item{locs}{a matrix of locations, each row specifying
  two locations to find the distance between. The first two
  columns are the "from" longitude and latitude,
  respectively; the third and fourth columns are the "from"
  longitude and latitude, respectively.}

  \item{dist}{logical: return the distances as a
  \code{\link{dist}} object (distance matrix)?}
}
\value{
either a matrix or a distance matrix containing the
distances between all the points.
}
\description{
Find the distance between locations specified by
latitude/longitude in degrees.
}

