pointOnLine <- function (cool, coop) {
    nearest_points = vapply(2:nrow(cool),
                            function(x) pointOnSegment(cool[(x - 1):x, ], coop),
                            FUN.VALUE = c(0, 0, 0, 0))
    dchain <- sqrt(diff(cool[,1])^2 + diff(cool[,2])^2)
    chainage <- c(0,cumsum(dchain))
    id <-  which.min(nearest_points[3, ])
    nearest_points[4,id] <- nearest_points[4,id] + chainage[id]
    nearest_points[, id]                    # 'x','y','d','chain'
} # end function pointOnLine
