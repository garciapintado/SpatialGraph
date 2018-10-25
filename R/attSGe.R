attSGe <- function(SG, att, eID, val, default) {
  SG@e@data[,att] <- default
  eids <- match(eID, row.names(SG@e@data))
  SG@e@data[eids,att] <- val
  return(SG)
} # end function attSGe
