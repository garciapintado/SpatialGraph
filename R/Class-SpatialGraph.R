setClass("SpatialGraph",
         representation(v='SpatialPointsDataFrame',
                        e='SpatialLinesDataFrame',
                        dist='matrix',
                        path='matrix'),
         validity = function(object) {
           if (!all(object@e@data[,'v0'] %in% object@v@data[,'ID']))
             stop("v0 vertices in edges not found")
           if (!all(object@e@data[,'v1'] %in% object@v@data[,'ID']))
             stop("v1 vertices in edges not found")
           return(TRUE)
         }
)
