#' This function takes input to the lagr or lagr.tune function, returning a model matrix, coordinates, weights, and response.
#' 
#' 
lagr.parse.model.frame = function(formula, data, weights, coords, fit.loc, longlat, na.action, contrasts) {
    # Determine whether coords was supplied as a character vector or as a symbolic expression
    coords.is.char = FALSE
    try(coords.is.char <- is.character(coords), silent=TRUE)
    
    # If coords was supplied as a character vector, then mave sure we include the named variables in the model.frame
    if (coords.is.char)
        vars = unlist(c("formula", "data", "weights", "na.action"), coords)
    else
        vars = c("formula", "data", "weights", "na.action")
    
    # Match the variables that are referenced in the function call
    mf = match.call(expand.dots=FALSE)
    m <- match(vars, names(mf), 0)
    mf <- mf[c(1, m)]
    
    # We need to evaluate the function call in the next environment up because we've added a level of abstraction with this function.
    mf[[1]] <- as.name("model.frame")
    mf <- eval.parent(substitute_q(mf, env=sys.frame(sys.parent())))
    mt <- attr(mf, "terms")

    # If the data was provided as a spatial data frame, then extract both the data and the coordinates.
    if (is(data, "Spatial")) {
        if (!missing(coords)) 
            warning("data is Spatial* object, ignoring coords argument")
        coords <- coordinates(data)
        if ((is.null(longlat) || !is.logical(longlat)) && !is.na(is.projected(data)) && !is.projected(data))
            longlat <- TRUE
        else longlat <- FALSE
    } else {
        # Make sure coordinates were specified
        if (missing(coords)) 
            stop("Observation coordinates have to be given")
        
        # Get the coords from the data:
        if (coords.is.char) coords = data[,coords]
        else {
            coords.expression = substitute(coords, env=sys.frame(sys.parent()))
            coords.expression[[1]] = as.name('cbind')
            coords = eval(coords.expression, data)
        }
        
        # Only interpret the coordinates as latitude/longitude values if the longlat variable is TRUE
        if (is.null(longlat) || !is.logical(longlat)) 
            longlat <- FALSE
    }
    
    # Get the matrices of distances and weights
    D.coords = rbind(coords, fit.loc)
    n = nrow(D.coords)
    if (longlat) {
        D = as.matrix(earth.dist(D.coords),n,n)
    } else {
        Xmat = matrix(rep(D.coords[,1], times=n), n, n)
        Ymat = matrix(rep(D.coords[,2], times=n), n, n)
        D = sqrt((Xmat-t(Xmat))**2 + (Ymat-t(Ymat))**2)
    }
    
    # Find the possible range of bandwidths (for use with the adaptive bandwith methods - knn or nen)
    bbox <- cbind(range(coords[, 1]), range(coords[, 2]))
    difmin <- spDistsN1(bbox, bbox[2, ], longlat)[1]
    if (any(!is.finite(difmin))) 
        difmin[which(!is.finite(difmin))] <- 0
    min.dist = difmin / 300
    max.dist = 10 * difmin
    
    # Get the data and the weights
    y <- model.response(mf, "numeric")
    x <- model.matrix(mt, mf, contrasts)
    w <- model.weights(mf)

    # Check for problems with the (prior) weights
    if (!is.null(w) && !is.numeric(w)) 
        stop("'weights' must be a numeric vector")
    if (is.null(w)) 
        w <- rep(1, nrow(data))
    if (any(is.na(w))) 
        stop("NAs in weights")
    if (any(w < 0)) 
        stop("negative weights")
    
    return(list(x=x, y=y, w=w, coords=coords, dist=D, max.dist=max.dist, min.dist=min.dist, mt=mt))
}