epanechnikov = function(R, bw) {
  ifelse( R < bw, 1-(R/bw)**2, 0)
}

cubic = function(d, bw) {
  ifelse(d<bw, 1 - (7*(d/bw)**2 - 8.75*(d/bw)**3 + 3.5*(d/bw)**5 - 0.75*(d/bw)**7), 0)
}

spherical = function(d, bw) {
  ifelse(d<bw, 1 - 1.5*(d/bw) + 0.5*(d/bw)**3, 0)
}

bisquare = function(R, bw) {
  ifelse( R < bw, (1 - (R/bw)**2)**2, 0)
}
