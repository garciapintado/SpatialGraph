polylineChainage <- function(xy) {
  # purpose
  #   obtain the chainage of nodes along a polyline [2-col matrix]
    dchain <- sqrt(diff(xy[,1])^2 + diff(xy[,2])^2)
    chainage <- c(0,cumsum(dchain))
}
