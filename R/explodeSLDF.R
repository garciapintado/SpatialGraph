explodeSLDF <- function(SLDF, FID = NULL) {
 # +++ purpose +++
 # a Lines slot in a SpatialLinesDataFrame may have several Line slots.
 # This is function is not fast, but it provides a way of obtaining a SpatialLinesDataFrame
 # in which each Lines contains just one Line. The process is done by creating a Lines slot for each Line, and duplicating
 # the attributes data.frame, with the exception of the unique identifier which increases with every newly created Lines

 # the feature identifiers for new Lines will be automatically set as a character representation of sequentially increasing numeric identifiers (matching the rownames of the attribute table).

 IDs <- sapply(SLDF@lines, function(x) slot(x, 'ID'))
 IDsNum <- as.numeric(IDs)
 idnum <- max(IDsNum, na.rm=TRUE)+1
 if (idnum < 0)                                             # all character IDs
   idnum <- 0

 if (!is.null(FID)) {
   FIDs <- SLDF@data[,FID]
   if (!is.numeric(FIDs))
     warning('explodeSLDF: FID is not a numeric field')
   FIDs   <- as.numeric(FIDs)
   fidnum <- max(FIDs, na.rm=TRUE)
   if (fidnum < 0)                                             # all character IDs
     fidnum <- 0
 }

 i <- 1
 while (i <= length(SLDF)) {
   nLines <- length(SLDF@lines[[i]]@Lines)
   if (nLines > 1) {
     SLDFsub <- SLDF[i,]
     SLDF@lines[[i]]@Lines <- SLDF@lines[[i]]@Lines[1]
     for (j in 2:nLines) {
       idnum  <- idnum + 1
       fidnum <- fidnum + 1
       SLDFtmp <- SLDFsub
       SLDFtmp@lines[[1]]@Lines <- SLDFtmp@lines[[1]]@Lines[j]
       SLDFtmp@lines[[1]]@ID  <- as.character(idnum)
       rownames(SLDFtmp@data) <- as.character(idnum)
       SLDFtmp@data[,FID] <- fidnum
       SLDF <- rbind(SLDF,SLDFtmp)
     }
   }
   i <- i + 1
 }
 return(SLDF)
} # end function explodeSLDF
