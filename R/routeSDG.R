routeSDG <- function(SDG, FUN='cumsum', ifld='inflow') {
 # SDG :: SpatialGraph [considreded as a directed graph] with slots
 #         @v         'SpatialPointsDataFrame'
 #         @e         'SpatialLinesDataFrame'

 # routeSDG assumes each vertex can provide an input to the system
 # the function returns a 'flow' vector.
 # in input, every vertex can be a source/sink, junction, or junction+source/sink
 # SDG[[1]] is the vertex information, and it must have the fields:
 #   coordinates
 #   inflow
 #   ID
 #  where
 #   inflow  > 0 in input vertex is a source, or junction+source
 #   inflow  < 0 in input vertex is a sink, or junction+sink
 #   inflow == 0 in input vertex is a simple junction
 #
 # SDG[[2]]@data provides the routing information, with at least the fields
 #  v0 : from vertex 'ID'
 #  v1 : to   vertex 'ID'
 #  div: diversion: 1 for single downstream propagation, < 1 for diversions, where the sum of downstream flows [normally]
 #       would sum to 1

 # FUN
 #  cumsum: input is instantly accumulated downstream [the only option currently available]
 #
 # require(sp)

 ldf <- as.data.frame(SDG@e)
 nv <- length(SDG@v)
 todo <- rep(TRUE,nv)
 flow <- rep(0,nv)

 iw <- 0
 while (TRUE) {
   iw <- iw + 1
   if (iw == 1) {
     ids <- !(SDG@v@data[,'ID'] %in% ldf[,'v1'])                                             # vertex iterative index
     flow[ids] <- SDG@v@data[ids,ifld]
     todo[ids] <- FALSE
   }

   for (i in 1:nv) {
     if (todo[i]) {
       ild  <- ldf[,'v1'] == SDG@v@data[i,'ID']
       div  <- ldf[ild,'div']
       iup <-  ldf[ild,'v0']                                                         # upstream vertices
       iup <- match(iup, SDG@v@data[,'ID'])
       if (all(!todo[iup])) {
         if (FUN == 'cumsum') {
           flow[i] <- sum(flow[iup] * div) + SDG@v@data[i,ifld]
           todo[i] <- FALSE
         }
       }
     }
   }
   if (all(!todo))
    break
 }
 return(flow)
} # end function routeSDG
