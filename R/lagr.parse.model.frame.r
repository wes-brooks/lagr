#' This function takes input to the lagr or lagr.tune function, returning a model matrix, coordinates, weights, and response.
#' 
#' 
lagr.parse.model.frame = function(formula, data, family, weights, coords, fit.loc, longlat, na.action, contrasts) {
    # Get the exponential family of the response
    if(is.character(family))
        family <- get(family, mode = "function", envir = parent.frame())
    if(is.function(family)) family <- family()
    if(is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
    }
    
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
        if (coords.is.char) coords = cbind(data[,coords])
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
    coords = as.matrix(coords)
    q = ncol(coords)
    D.coords = rbind(coords, fit.loc)
    n = nrow(D.coords)

    #Two dimensional effect modifying parameter:
    if (q==2 && longlat) {
        #If data was specified in terms of latitude/longitude, use earth distance
        D = as.matrix(earth.dist(D.coords),n,n)
    } else {
        #Otherwise, use the pythagorean distance
        squaredist = matrix(0,n,n)
        D.coords = as.matrix(D.coords)
        
        for (c in 1:ncol(as.matrix(coords))) {
            this.coord = matrix(rep(D.coords[,c], times=n), n, n)
            squaredist = squaredist + (this.coord-t(this.coord))**2
        }
        D = sqrt(squaredist)
    }
    
    max.dist = 10 * max(D)
    min.dist = max.dist / 3000
    
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
    
    return(list(x=x, y=y, w=w, family=family, coords=coords, dist=D, max.dist=max.dist, min.dist=min.dist, mt=mt, dim=q))
}