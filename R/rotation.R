rotation <- function (coords, radian) {
  # +++ purpose +++
  # rotate points counter clockwise
  # coords: 2-col matrix of [x,y] coordinates
  # radian: angle to conduct counterclokwise rotation (positive angles)

  x <- c(cos(radian), -sin(radian))
  y <- c(sin(radian), cos(radian))
  nlecoord = coords %*% cbind(x, y)
  return(nlecoord)
} # end function rotation
