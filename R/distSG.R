distSG <- function(SG, x, y = NULL, euc = TRUE, wei = NULL, getpath = FALSE) {
 # calculate along-network distance between sparse points
 # SG :: SpatialGraph object
 # x  :: from-point SpatialPointsDataFrame with the fields 'chain' and 'eID', providing its basic relationship with the SpatialGraph
 # y  :: to-point SpatialPointsDataFrame with the fields 'chain' and 'eID', providing its basic relationship with the SpatialGraph
 # euc:: whether to consider the Euclidean distance as minimum threshold for along-network distance
 # wei:: rules for weighting
 #
 # The distance between 2 points is defined within this function
 # as the minimum along-network distance. The distance itself between the points
 # in x,y and the network is neglected in the function. But the maximum of both
 # the Euclidean distance and the along-network distance is returned if euc=TRUE [default]
 #require(pracma)

  ldf    <- SG@e@data
  eID    <- sapply(SG@e@lines, function(x) slot(x, 'ID'))
  x$v0   <- ldf$v0[match(x$eID,eID)]        # match(x$eID,eID) maps edge identifier into edge iterative order
  x$v1   <- ldf$v1[match(x$eID,eID)]
  x$elen <- ldf$len[match(x$eID,eID)]

  if (is.null(y)) {
    y <- x
  } else {
    y$v0   <- ldf$v0[match(y$eID,eID)]
    y$v1   <- ldf$v1[match(y$eID,eID)]
    y$elen <- ldf$len[match(y$eID,eID)]
  }
  # along network distance matrix
  nx <- nrow(x)
  ny <- nrow(y)

  rid <- matrix(NA,nx,4)
  cid <- matrix(NA,ny,4)

  disa <- array(NA,dim=c(nx,ny,5))
  rid[,1] <- match(x$v0, rownames(SG@dist))                             # rows within the vertex distance matrix
  cid[,1] <- match(y$v0, rownames(SG@dist))                             # columns "         "        "        "
  disa[,,1] <- SG@dist[rid[,1],cid[,1]] + rep(x$chain,ny) + rep(y$chain, each=nx)
  rid[,2] <- match(x$v0, rownames(SG@dist))                             #  "
  cid[,2] <- match(y$v1, rownames(SG@dist))                             #  "
  disa[,,2] <- SG@dist[rid[,2],cid[,2]] + rep(x$chain,ny) + rep(y$elen - y$chain, each=nx)
  rid[,3] <- match(x$v1, rownames(SG@dist))                             #  "
  cid[,3] <- match(y$v0, rownames(SG@dist))                             #  "
  disa[,,3] <- SG@dist[rid[,3],cid[,3]] + rep(x$elen - x$chain,ny) + rep(y$chain, each=nx)
  rid[,4] <- match(x$v1, rownames(SG@dist))                             #  "
  cid[,4] <- match(y$v1, rownames(SG@dist))                             #  "
  disa[,,4] <- SG@dist[rid[,4],cid[,4]] + rep(x$elen - x$chain,ny) + rep(y$elen - y$chain, each=nx)

  matche <- matrix(x$eID,nx,ny) == matrix(y$eID,nx,ny,byrow=TRUE)
  edgs <- matrix(x$eID,nx,ny)[matche]
  disa[,,5] <- abs(rep(x$chain,ny) - rep(y$chain, each=nx))
  disa[,,5][!matche] <- Inf
  dis <- apply(disa,c(1,2),min)

  if (euc) {
    dise <- distmat(coordinates(x),coordinates(y))
    dis <- pmax(dis,dise)
  }
  if (getpath)
   pathmat <- matrix(list(),nx,ny)
     
  if (!is.null(wei)) {   # also obtain state-dependent weights and return results list
    if (nrow(SG@path) != length(SG@v))
      stop('distSG -- ERR01 - no slot of name "path" for this "SpatialGraph" object')
    diw <- apply(disa,c(1,2),which.min)               # which min path out of the 5 options
    calcw <- diw != 5 & !is.infinite(dis) & dis > 0.0 # which paths will add weight
    pathid <- array(NA,dim=c(nx,ny,2))
    wemat  <- matrix(1,nx,ny)
    for (i in 1:nx) {
      for (j in 1:ny) {
        if (calcw[i,j]) {
          # reverse <- FALSE
          pathid[i,j,1] <- rid[i,diw[i,j]] # row indexes within SG@dist
          pathid[i,j,2] <- cid[j,diw[i,j]] # col indexes within SG@dist
          # if (pathid[i,j,1] > pathid[i,j,2]) { # path is upper triangular matrix
          #   pathid[i,j,] <- pathid[i,j,2:1]
          #   reverse <- TRUE
          # }
          eIDs <-   SG@path[[pathid[i,j,1],pathid[i,j,2]]]$e
          # if (reverse)
          #   eIDs <- rev(eIDs)
          eIDs <- c(x$eID[i],eIDs,y$eID[j])
          westa <-  ldf[,wei][match(eIDs,eID)]
          westa <- westa[-length(westa)]/westa[-1]
          westa[westa < 1.] <- 1./westa[westa<1]
          westa[westa>10.] <- 10. # hard-coded threshold for single state-dependent weights in the chain
          wemat[i,j] <- prod(westa)
          if (getpath)
            pathmat[[i,j]] <- eIDs 
        } # if calcw
      } # for j
    } # for i
    if (!getpath)
      dis <- list(dis=dis,wei=wemat)
    else
      dis <- list(dis=dis,wei=wemat,eID=pathmat)
  }
  return(dis)
} # end function distSG
