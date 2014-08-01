\name{deg.dist}
\alias{deg.dist}
\title{Find the distance from the "from" location to a vector of "to" locations, given their latitude and longitude. Used by the earth.dist function.}
\usage{
deg.dist(long1, lat1, long2, lat2)
}
\arguments{
  \item{long1}{the "from" longitude, in degrees}

  \item{lat1}{the "from" latitude, in degrees}

  \item{long2}{the vector of "to" longitudes, in degrees}

  \item{lat2}{the vector of "to latitudes, in degrees}
}
\value{
numeric scalar distance from the "from" location to the
"to" locations.
}
\description{
Find the distance from the "from" location to a vector of
"to" locations, given their latitude and longitude. Used by
the earth.dist function.
}

