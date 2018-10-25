SpatialGraph <- function(v, e, dist=NULL, path=NULL) {
  if (is.null(dist))
    dist <- as.matrix(NA)
  if (is.null(path))
    path <- matrix(list())
  new("SpatialGraph", v=v, e=e, dist=dist, path=path)
}


