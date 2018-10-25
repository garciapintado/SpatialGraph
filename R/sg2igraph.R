sg2igraph <- function(sg, directed=FALSE) {
  # map a SpatialGraph object into an igraph object
  # a SpatialGraph is a S4 class object with the slots:
  #  @v    : vertices [SpatialPointsDataFrame]
  #  @e    : edges    [SpatialLinesDataFrame]
  #  @dist : across-network distance matrix among vertices
  #require(igraph)
  g <- graph.empty(directed=directed)
  nv <- nrow(sg@v)
  spID <- sg@v@data[,'ID']
  g <- add.vertices(g, nv, x=coordinates(sg@v)[,1], y=coordinates(sg@v)[,2], spID=spID)
  #iID <- as.numeric(V(g))
  ldf <- as.data.frame(sg@e)
  ifrom  <- match(ldf[,'v0'],spID)
  ito    <- match(ldf[,'v1'],spID)
  iedges <- matrix(c(ifrom,ito), ncol=2)
  g <- add.edges(g, t(iedges), spfrom=ldf[,'v0'], spto=ldf[,'v1'],
                 div=ldf[,'div'], len=ldf[,'len'])
  return(g)
} # end function sg2igraph
