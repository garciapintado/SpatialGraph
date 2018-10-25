splitSLDF <- function(SLDF, SPDF, dmax=NULL) {
 # split 1-Line Lines in a SpatialLinesDataFrame, by intersection with a set of points
 #

 if (class(SLDF) == 'SpatialLines')
   cls <- 'sl'
 else if (class(SLDF) == 'SpatialLinesDataFrame')
   cls <- 'sldf'
 else
   stop ('splitSLDF - ERR01 - check input object class')

 xyp <- coordinates(SPDF)
 if (is.null(dmax))
   dmax <- mean(diff(apply(xyp,2,range))) / 1000  # e.g.== 10m / 10km

 ne <- length(SLDF)
 acrn <- 0                             # accumulated row number

 for (ie in 1:ne) {
   xy  <- coordinates(SLDF[ie,])[[1]]
   if (length(xy) > 1)
    stop('splitSLDF: just accept 1-Line Lines in input SpatialLinesDataFrame objects')
   xy <- xy[[1]]
   xy  <- splitPolyline(xy, xyp, dmax=dmax) # list with individual transects as polylines
   nLines <- length(xy)
   e <- vector('list', nLines)
   if (cls == 'sldf')
     edf <- SLDF@data[ie,]
   for (i in 1:nLines) {
     acrn <- acrn + 1
     e[[i]] <- Line(xy[[i]])
     e[[i]] <- Lines(e[i], ID=acrn)
     if (i > 1 && cls == 'sldf')
       edf <- rbind(edf, SLDF@data[ie,])
   }

   e <- SpatialLines(e)
   if (cls == 'sldf') {
     row.names(edf) <- sapply(e, function(x) slot(x, 'ID'))
     e <- SpatialLinesDataFrame(e, edf)
   }
   if (ie == 1)
     SLDFou <- e
   else
     SLDFou <- rbind(SLDFou,e)
 }
  return(SLDFou)
} # end function splitSLDF
