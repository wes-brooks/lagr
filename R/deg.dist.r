#' Find the distance from the "from" location to a vector of "to" locations, given their latitude and longitude. Used by the earth.dist function.
#' 
#' @param long1 the "from" longitude, in degrees
#' @param lat1 the "from" latitude, in degrees
#' @param long2 the vector of "to" longitudes, in degrees
#' @param lat2 the vector of "to latitudes, in degrees
#' 
#' @return numeric scalar distance from the "from" location to the "to" locations.
deg.dist = function (long1, lat1, long2, lat2) {
    rad <- pi/180
    a1 <- lat1 * rad
    a2 <- long1 * rad
    b1 <- lat2 * rad
    b2 <- long2 * rad
    dlon <- b2 - a2
    dlat <- b1 - a1
    a <- (sin(dlat/2))^2 + cos(a1) * cos(b1) * (sin(dlon/2))^2
    c <- 2 * atan2(sqrt(a), sqrt(1 - a))
    R <- 40041.47/(2 * pi)
    d <- R * c
    return(d)
}