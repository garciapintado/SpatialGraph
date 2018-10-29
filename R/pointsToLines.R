pointsToLines <- function (points, lines, withAttrs = TRUE, withDis = TRUE, withChain = TRUE) {
    # browser()
    #require("rgeos")
    # points :: SPDF, SP, or 2-col matrix
    # lines  :: SpatialLinesDataFrame
    if (!is(points, "SpatialPointsDataFrame")) { # i.e. 'SpatialPoints' or 'matrix'
      if (missing(withAttrs))
        withAttrs = FALSE
      else if (withAttrs == TRUE)
        stop("withAttrs = TRUE just available for SpatialPointsDataFrame objects")
    }
    #if (!is.na(maxDist)) {
    #    w = gWithinDistance(points, lines, dist = maxDist, byid = TRUE)
    #    validPoints = apply(w, 2, any)
    #    validLines = apply(w, 1, any)
    #    points = points[validPoints, ]
    #    lines = lines[validLines, ]
    #}
    if (is(points,'matrix')) {
      coordsPoints <- points[,1:2,drop=FALSE]
      points <- as.data.frame(coordsPoints)
      names(points) <- c('x','y')
      coordinates(points) <- c('x','y')
    } else {
      coordsPoints = coordinates(points)
    }
    n = length(points)
    d = gDistance(points, lines, byid = TRUE)                                          # [m,n] matrix, m= number of lines
    nearest_line_index = apply(d, 2, which.min)
    d = apply(d, 2, min)

    coordsLines = coordinates(lines)

    mNewCoords = vapply(1:n,                                                                          # [4,n] matrix : x,y [crossing point], d [distance point-crossing point], chain [differential chainage]
                        function(x) pointOnLine(coordsLines[[nearest_line_index[x]]][[1]],
                                                coordsPoints[x,,drop=FALSE]),
                        FUN.VALUE = c(0, 0, 0, 0))                                                    #
    #if (!is.na(maxDist))
    #  nearest_line_id = as.numeric(rownames(d)[nearest_line_index]) + 1
    #else
    #  nearest_line_id = nearest_line_index
    eID <- sapply(lines@lines, function(x) slot(x, 'ID'))
    if (withAttrs)
      df = cbind(points@data, lid=nearest_line_index, eID = eID[nearest_line_index], stringsAsFactors=FALSE)
    else
      df = data.frame(lid=nearest_line_index, eID = eID[nearest_line_index], stringsAsFactors=FALSE)
    if (withDis)
      df = cbind(df, dis=d)
    if (withChain)
      df = cbind(df, chain=mNewCoords[4,])

    SpatialPointsDataFrame(coords = t(mNewCoords[1:2,]), data = df,
                           proj4string = CRS(proj4string(lines)))

} # end function pointsToLines
