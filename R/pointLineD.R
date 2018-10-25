pointLineD <- function(xy, xyp) {
 # +++ purpose +++
 # obtain the closest point [to] and distance [from] a line segment to a set of point/s
 # xy  (I) : 2 x 2 matrix of coordinates defining the line segment
 # xyp (I) : p x 2 matrix of the external points to estimate distances to the line

 # return a list with the components:
 # xy : matrix [2 x 2], with the xy components of the
 #      segment starting point
 #      segment end point
 #      given point [xyp] to obtain distance
 #      perperdicular projection of xyp onto the line [may be out of the line segment]
 # d  : point-line distance
 # dl : distance, along the segment, between the first point in segment and the projection of xyp onto the line
 # cro: boolean [coded as 0s/1s] indicating whether the projected point [xyc] falls within the segment
    
 alpha <- atan(abs((xy[2,2]-xy[1,2])/(xy[2,1]-xy[1,1])))
 if (xy[1,1] >  xy[2,1] && xy[1,2] <= xy[2,2])                       # 2nd cuadrant
   alpha <- pi - alpha
 if (xy[1,1] >  xy[2,1] && xy[1,2] > xy[2,2])                        # 3rd cuadrant
   alpha <- alpha + pi
 if (xy[1,1] <= xy[2,1] && xy[1,2] >= xy[2,2])                       # 4th cuadrant
   alpha <- 2*pi - alpha

 if (length(xyp) == 2)
   xyp <- matrix(xyp,ncol=2)
 p <- nrow(xyp)

 # move into the rotated space
 xyr <- rotation(rbind(xy,xyp), pi/2-alpha)                 # counter clockwise rotation of segment & source points
 xyc <- matrix(c(rep(xyr[1,1],p),xyr[3:(2+p),2]), ncol=2)  # rotated coordinates of the line points nearest to xyp points (normal projections)
 x0 <- xyr[1,1]
 y0 <- xyr[1,2]
 y1 <- xyr[2,2]
 d     <- abs(x0 - xyr[3:(2+p),1])          # points-line distance
 dc    <- xyc[,2] - y0                      # diferential chainage over [x0,y0] (> 0 if the projection goes in the segment direction)
 cross <- rep(0,p)
 cross[xyc[,2] >= y0 & xyc[,2] <= y1] <- 1  # 1 = closest to line than to nodes

 # rotate back into physical space
 xy  <- rotation(xyr, alpha-pi/2)
 xyc <- rotation(xyc, alpha-pi/2)
 ans <- list()
 ans$xy    <- round(xy[1:2,],10)
 ans$xyp   <- round(xy[3:(2+p),,drop=FALSE],10)
 ans$xyc   <- xyc
 ans$d     <- d
 ans$dc    <- dc
 ans$cross <- cross
 return(ans)
} # end function pointLineD
