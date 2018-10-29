pointOnSegment <- function (s, p) {
    # s :: [2,2] matrix, each row representing one of the end points in the segment
    # p :: [2]   vector, representing the X,Y coordinates of the point

    ap = c(p[1] - s[1, 1], p[2] - s[1, 2])
    ab = c(s[2, 1] - s[1, 1], s[2, 2] - s[1, 2])            # [dx,dy]
    t = sum(ap * ab)/sum(ab * ab)
    t  = ifelse(t < 0, 0, ifelse(t > 1, 1, t))
    x = s[1, 1] + ab[1] * t
    y = s[1, 2] + ab[2] * t
    dc = ifelse(t == 0, 0, ifelse(t == 1, sqrt(sum(ab^2)),
                                          sqrt((x - s[1,1])^2 + (y - s[1,2])^2)))
    result = c(x, y, sqrt((x - p[1])^2 + (y - p[2])^2), dc)
    names(result) = c("x", "y", "d", "chain")
    result
} # end function pointOnSegment
