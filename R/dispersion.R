dispersion <- function(y, mu, w, df, family) {
    sumw = sum(w)
    
    dev.resid.values = family$dev.resids(y, mu, w)
    
    if (family$family %in% c("poisson", "binomial")) {dispersion = 1}
    else if (sumw > 0) {dispersion = sum(resid) / (sumw-df)}
    else {dispersion = 1}
    
    return(dispersion)
}