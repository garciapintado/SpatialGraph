splitPolyline <- function(xy, xyp, dmax) {
  # xy:  polyline coordinates [2-col matrix]
  # xyp: points to split the polyline [2-col matrix]

  # disregard points representing the polyline limits
  ncoo <- nrow(xy)
  np   <- nrow(xyp)

  ver0 <-  rowSums(abs(xyp - matrix(xy[1,],byrow=TRUE,ncol=2,nrow=np))) == 0
  ver1 <-  rowSums(abs(xyp - matrix(xy[ncoo,],byrow=TRUE,ncol=2,nrow=np))) == 0
  isver <- ver0 | ver1
  ver0 <-  xyp[ver0,,drop=FALSE]
  ver1 <-  xyp[ver1,,drop=FALSE]
  xyp <- xyp[!isver,,drop=FALSE]

  if (length(xyp) == 0) {
    xyl <- list(xy)
    return(xyl)
  }
  # obtain crossing points and clip new polylines to the input points
  xyp <- unique(xyp)                              # remove possible duplicates
  pcha <- polylineChainage(xy)                    # polyline chainage at nodes
  plen <- polylineLength(xy)                      # polyline length
  xyc  <- pointPolylineD(xy,xyp)                  # chainage table

  isfar <- xyc[,'dis'] > dmax
  xyc  <- xyc[,'chain0'] + xyc[,'dc']             # chainages
  xyp <- xyp[!isfar,,drop=FALSE]
  xyc <- xyc[!isfar]
  if (length(xyc) == 0) {
    xyl <- list(xy)
    return(xyl)
  }

  isou <- xyc == 0 | xyc >= plen
  xyp <- xyp[!isou,,drop=FALSE]
  xyc <- xyc[!isou]
  if (length(xyc) == 0) {
    xyl <- list(xy)
    return(xyl)
  }
  isdu <- duplicated(xyc)
  xyp <- xyp[!isdu,,drop=FALSE]
  xyc <- xyc[!isdu]

  ord <- order(xyc)
  xyp <- xyp[ord,,drop=FALSE]
  xyc <- xyc[ord]

  nc <- length(xyc)
  xyl <- vector('list',nc+1)
  interv <- findInterval(pcha,xyc) + 1
  nt <-  length(xyc)+1                                     # number of transects
  for (i in 1:nt) {
    ids <- interv == i
    xyl[[i]] <- xy[ids,,drop=FALSE]
    if (i > 1)                                    # clip
      xyl[[i]] <- rbind(xyp[i-1,],xyl[[i]])
    if (i < nt)
      xyl[[i]] <- rbind(xyl[[i]],xyp[i,])
    if (i == 1)
      xyl[[i]] <- rbind(ver0,xyl[[i]])
    if (i == nt)
      xyl[[i]] <- rbind(xyl[[i]],ver1)
    isdu <- duplicated(xyl[[i]])
    xyl[[i]] <- xyl[[i]][!isdu,,drop=FALSE]
  }
  return(xyl)
} # end function splitPolyline
