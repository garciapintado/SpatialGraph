sgChVIDs <- function(obj, IDa, IDp = NULL) {
  nv <- nrow(slot(obj, "v"))
  if (length(IDa) != nv)
    stop("lengths differ")
  if (length(IDa) > length(unique(IDa)))
    stop("duplicate IDa")
  if (is.null(IDp))
    IDp <- obj@v@data[,'ID']

  vs <- obj@v@data[,'ID']
  id <- match(vs,IDp)
  obj@v@data[,'ID'] <- IDa[id]

  ev <- c('v0','v1')
  for (i in 1:2) {
    vs <- obj@e@data[,ev[i]]
    id <- match(vs,IDp)
    obj@e@data[,ev[i]] <- IDa[id]
  }
  return(obj)
} # end function sgChVIDs
