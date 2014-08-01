#' Compute \code{y*log(y)}
ylogy = function(y) {
    return( ifelse(y==0, 0, y*log(y)) )
}