\name{pointOnSegment}
\alias{pointOnSegment}
\title{Snap a points to a segment}
\description{
  This function snaps a point to a segment based on the
  minimum distance between the point and the segment}
\usage{
pointOnSegment(s, p)
}
\arguments{
  \item{s}{[2,2] matrix giving the coordinates of the line, one point
    per row}
  \item{p}{2-length vector repsenting the point}
}
\value{
 A 4-length vector, with 'x','y' [coordinates of the point snapped to
 the segment], 'd' [distance from the input point to the new snapped
 point], and 'chain' [distance from the starting
 of the segment to the snapped point]
}
\author{Javier Garcia-Pintado}
\seealso{\code{\link[sp]{Spatial-class}}}
\keyword{ utilities }
\keyword{ spatial }
