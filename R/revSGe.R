revSGe <- function(SG, eID) {
 # reverse the coordinates of the polylines in a SpatialGraph
 FID <- row.names(SG@e@data)
 lid <- match(eID,FID)
 for (i in lid) {
     coo <-  coordinates(SG@e[i,])[[1]][[1]]
     coo <- coo[nrow(coo):1,]
     SG@e@lines[[i]]@Lines <- list(Line(coo))
     SG@e@data[i,c('v0','v1')] <- SG@e@data[i,c('v1','v0')]
  }
  return(SG)
} # end function revSGe
