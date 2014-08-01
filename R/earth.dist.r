#' Find the distance between locations specified by latitude/longitude in degrees.
#' 
#' @param locs a matrix of locations, each row specifying two locations to find the distance between. The first two columns are the "from" longitude and latitude, respectively; the third and fourth columns are the "from" longitude and latitude, respectively.
#' @param dist logical: return the distances as a \code{\link{dist}} object (distance matrix)?
#' 
#' @return either a matrix or a distance matrix containing the distances between all the points.
earth.dist = function (locs, dist=TRUE) 
{
    if (class(locs) == "SpatialPoints") 
        locs <- coordinates(locs)
    
    name <- list(rownames(locs), rownames(locs))
    n <- nrow(locs)
    z <- matrix(0, n, n, dimnames = name)
    z <- outer(1:n, 1:n, function(i, j) deg.dist(long1=locs[i,1], lat1=locs[i, 2], long2=locs[j, 1], lat2=locs[j,2]))

    if (dist == TRUE) 
        z <- as.dist(z)
    
    return(z)
}