textSGe <- function(SG, acol='wheat', tcol='navyblue', arr.length=0.4) {

  polylineChainage <- function(xy) {
    # +++ purpose +++
    #   obtain the chainage of nodes along a polyline [2-col matrix]
    dchain <- sqrt(diff(xy[,1])^2 + diff(xy[,2])^2)
    chainage <- c(0,cumsum(dchain))
  }
  
  ne <- length(SG@e)
  einf <- vector('list',ne)
  ldf <- SG@e@data
  IDs <- row.names(ldf)
  for (i in 1:ne) {
    coo  <- coordinates(SG@e[i,])[[1]][[1]]
    chai <- polylineChainage(coo) 
    midc <- max(chai)/2
    lid  <- max(which(chai <= midc))
    hid  <- which(chai > midc)[1]
    coo  <- coo[c(lid,hid),]
    chai <- chai[c(lid,hid)]
    mcoo <- c(approx(x=chai,y=coo[,1],xout=midc)$y,     # x
              approx(x=chai,y=coo[,2],xout=midc)$y)     # y
    dcoo <- diff(coo)
    angl <- atan(dcoo[2]/dcoo[1]) * 180 / pi
    if (angl > 0 && dcoo[1] < 0)
      angl <- 180 + angl
    if (angl < 0) {
      if (dcoo[1] < 0)
        angl <- 180 + angl
      else
        angl <- 360 + angl
    }
    Arrowhead(mcoo[1],mcoo[2],angle=angl, lcol=acol, arr.col=acol, arr.type='triangle', arr.length=arr.length)
    text(mcoo[1],mcoo[2],label=IDs[i], pos=1,col=tcol)          
  }
  return(0)
} # end function textSGe
