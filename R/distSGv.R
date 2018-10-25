distSGv <- function(SG, getpath = FALSE) {
  iSG <- sg2igraph(SG)
  idist <- shortest.paths(iSG, weights=E(iSG)$len)             # along-network distance matrix
  colnames(idist) <- SG@v$ID
  rownames(idist) <- SG@v$ID
  SG@dist <- idist

  if (getpath) {
    nv <- length(SG@v)
    SG@path <- matrix(list(),nv,nv)
    for (iv in 1:nv) {
      toid <- rep(TRUE,nv)
      toid [is.infinite(SG@dist[iv,])] <- FALSE
      toid[0:(iv-1)] <- FALSE  # upper triangular matrix
      cc <- get.shortest.paths(iSG, from=iv, to=(1:nv)[toid], weights=E(iSG)$len, output='both')
      iv2i <- 0
      for (iv2 in (1:nv)[toid]) {
        iv2i <- iv2i + 1
        SG@path[[iv,iv2]] <- list(v=cc$vpath[[iv2i]],e=cc$epath[[iv2i]])
      }
    }
    uppath <- t(SG@path)[lower.tri(SG@path)]   # make symmetric [more storage but speeds up later computations]
    revlst <- function(x) {x <- lapply(x,FUN=rev)
                           if(length(x)==0)
                             x <- NULL
                           return(x)}
    uppath <- lapply(uppath,FUN=revlst)
    SG@path[lower.tri(SG@path)] <- uppath
  } # end if getpath
  return(SG)
} # end function distSGv
