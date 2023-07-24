sl2sg <- function(SL, clipd = NULL, getdist = TRUE, getpath = FALSE) {
  # purpose:
  #   map a network provided as SpatialLinesDataFrame into
  #   a SpatialGraph [SG]
  #
  # Details:
  #   a SG is a S4 class object with slots
  #   @v : SpatialPointsDataFrame representing the graph vertices
  #   @e : SpatialLinesDataFrame representing the graph edges
  #   @dist: across-network distance matrix between vertices
  #   @path: list with optimal-path matching the dist slot
  #   A difference with a standard [igraph] graph is that the edges may be polylines,
  #   where each polyline is represented by a SLDF@lines object, and its
  #   corresponding [1-row] entry in the SLDF@data data.frame
  #   The topology [connectivity] information is provided by SLDF@data
  #
  #   clipd: clipping distance: polylines within this distance of new junctions will be clipped by adding additional vertices
  #          TODO: add the clip to existing non-consistent junctions [very close start / end of lines]

  ## upgrade any individual Line to a lines slot

  if (!inherits(SL,'SpatialLines'))
    stop('sl2sg - ERR01 -- input object not SpatialLines nor SpatilLinesDataFrame')

  SL <- explodeSLDF(SL)
  ne <- length(SL)
  SL <- spChFIDs(SL,as.character(1:ne))

  ## get all polyline limits and intersections as vertices
  v <- NULL
  for (i in 1:length(SL)) {
    coo <-  coordinates(SL[i,])[[1]][[1]]
    v <- rbind(v, coo[c(1,nrow(coo)),])
  }
  for (i in 1:length(SL)) {
    for (j in 1:length(SL)) {
      if (i != j) {
        # coo <- gIntersection(SL[i,],SL[j,]) # rgeos::
        coo <- st_intersection(as(SL[i,], "sf"), as(SL[j,], "sf"))
        if (!is.null(coo))
          # v <- rbind(v,coordinates(coo))
          v <- rbind(v,st_coordinates(coo)) # sf::
      }
    }
  }
  v <- unique(v)
  v <- data.frame(v, row.names=NULL)
  names(v) <- c('x','y')
  v$inflow <- 0.0
  v$ID     <- 1:nrow(v)

  coordinates(v) <- c('x','y')

  nv <- nrow(v)

  if (is.null(clipd))
    clipd <- mean(diff(apply(coordinates(v),2,range))) / 1.0E04  # e.g.== 1m / 10km

  ## break polylines at junctions
  v_sf <- as(v, "sf")                             # "sf" "data.frame"
  deli <- rep(FALSE,ne)
  for (i in 1:ne) {
    # cat("i:", i, "\n")
    # get vertices intersecting the line
    xvid <- rep(FALSE, nv)
    #for (iv in 1:nv) {
      #xvid[iv] <- gWithinDistance(SL[i,], v[iv,], dist=clipd)   # rgeos:: crossing vertex ids
    #  
    close_ids <- st_is_within_distance(as(SL[i,], "sf"), v_sf, dist=clipd)[[1]]   # sf:: crossing vertex ids
    if (length(close_ids) > 0) {
      xvid[close_ids] <- TRUE   
    }
    if (sum(xvid) > 2) {
      deli[i] <- TRUE                             # mark for deletion
      SLpad <- splitSLDF(SL[i,],v[xvid,])
      FIDs <- length(SL) + 1:length(SLpad)
      SLpad <- spChFIDs(SLpad,as.character(FIDs))
      SL <- rbind(SL,SLpad)
    }
  }

  del <- rep(FALSE,length(SL))
  del[1:ne] <- deli
  SL <- SL[!del,]
  ne <- length(SL)
  SL <- spChFIDs(SL,as.character(1:ne))

  ## build topology
  v0  <- rep(0,ne)
  v1  <- rep(0,ne)
  len <- rep(0,ne)
  vcoo <- coordinates(v)
  for (i in 1:ne) {
    coo <- coordinates(SL[i,])[[1]][[1]]
    ncoo <- nrow(coo)
    lv   <- coo[c(1,ncoo),]

    v0[i] <- which(rowSums(abs(vcoo - matrix(lv[1,],byrow=TRUE,ncol=2,nrow=nv))) == 0)
    v1[i] <- which(rowSums(abs(vcoo - matrix(lv[2,],byrow=TRUE,ncol=2,nrow=nv))) == 0)
    len[i] <- polylineLength(coo)
  }
  newdat <- data.frame('v0'=v0,'v1'=v1,'div'=rep(1,ne),'len'=len)

  if (class(SL) == 'SpatialLinesDataFrame')
    SL@data <- cbind(SL@data,newdat)
  else # SpatialLines
    SL <- SpatialLinesDataFrame(SL,newdat)

  SG <- SpatialGraph(v=v, e=SL)

  iSG <- sg2igraph(SG)

  if (getdist) {
    idist <- shortest.paths(iSG, weights=E(iSG)$len)             # igraph:: along-network distance matrix
    colnames(idist) <- SG@v$ID
    rownames(idist) <- SG@v$ID
    SG@dist <- idist
  }
  if (getdist && getpath) {
    nv <- length(SG@v)
    SG@path <- matrix(list(),nv,nv)
    for (iv in 1:nv) {
      toid <- rep(TRUE,nv)
      toid [is.infinite(SG@dist[iv,])] <- FALSE
      toid[0:(iv-1)] <- FALSE
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
} # end function sl2sg
