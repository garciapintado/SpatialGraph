pointsPolylineD <- function(xy, xyp) {
  # +++ purpose +++
  # obtain the closest points in a polyline to a set of point/s
  # xy  (I) : n x 2 matrix of the segment coordinates which define the polyline
  # xyp (I) : p x 2 matrix of the external points to estimate distances to the line
  #
  # details: first the distance from the set of points to the lines defined by every single segment in the polyline
  # is obtained, then the distance to every single node in the polyline are also obtained. The lowest distance is chosen.

  if (class(xyp) == 'data.frame')
   xyp <- as.matrix(xyp)
  dimnames(xyp) <- NULL

  if (length(xyp) == 2)
    xyp <- matrix(xyp,ncol=2)
  p <- nrow(xyp)

  ncoo  <- nrow(xy)
  dchain <- sqrt(diff(xy[,1])^2 + diff(xy[,2])^2)
  chainage <- c(0,cumsum(dchain))

  crlst <- list() # list with point-lines crossing information
  d     <- matrix(NA,nrow=ncoo-1,p)
  cross <- matrix(NA,nrow=ncoo-1,p)
  for (i in 1:(ncoo-1)) {
    crlst[[i]] <- pointLineD(xy[c(i,i+1),],xyp)
    d[i,]     <- crlst[[i]]$d
    cross[i,] <- crlst[[i]]$cross
  }

  dord   <- apply(d, MARGIN=2, FUN=order)
  if (!is.matrix(dord))
    dord <- matrix(dord,ncol=p)
  crsort <- cross
  inode <- rep(0,p)
  for (ip in 1:p) {
    crsort[,ip] <- cross[dord[,ip],ip]
    inode[ip]   <- dord[which(crsort[,ip] == 1)[1],ip] # if no crossing normal segment is available
  }
  ansn <- c('inode','x0','y0','chain0','xc','yc','dc','dis')
  ans  <- as.data.frame(matrix(0, nrow=p, ncol=8, dimnames=list(NULL,ansn)))
  ans$inode  <- inode
  ans$x0     <- xy[inode,1]
  ans$y0     <- xy[inode,2]
  ans$chain0 <- chainage[inode]
  for (ip in 1:p) {
    ans$xc[ip]  <- c(crlst[[inode[ip]]]$xyc[ip,1],NA)[1]
    ans$yc[ip]  <- c(crlst[[inode[ip]]]$xyc[ip,2],NA)[1]
    ans$dc[ip]  <- c(crlst[[inode[ip]]]$dc[ip],   NA)[1]
    ans$dis[ip] <- c(crlst[[inode[ip]]]$d[ip],    NA)[1]
  }
  dis2nodes <- n2dist(xy,xyp)
  nodeclo   <- ans$dis > dis2nodes$dists # node closest than any segment
  nodeclo[is.na(nodeclo)] <- TRUE
  ans[nodeclo,'inode'] <- dis2nodes$neighs[nodeclo]
  ans[nodeclo,'dis']   <- dis2nodes$dist[nodeclo]
  ans[nodeclo,'x0'] <- ans[nodeclo,'xc'] <- xy[ans[nodeclo,'inode'],1]
  ans[nodeclo,'y0'] <- ans[nodeclo,'yc'] <- xy[ans[nodeclo,'inode'],2]
  ans[nodeclo,'chain0'] <- chainage[ans[nodeclo,'inode']]
  ans[nodeclo,'dc'] <- 0.0

  return(ans)
}
