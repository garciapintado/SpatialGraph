\name{pointOnLine}
\alias{pointOnLine}
\title{Snap a points to a line}
\description{
  This function snaps a point to a line based on the
  minimum distance between the point and the line}
\usage{
pointOnLine(cool, coop)
}
\arguments{
  \item{cool}{2-col matrix giving the coordinates of the line}
  \item{coop}{2-length vector repsenting the point}
}
\value{
 A 4-length vector, with 'x','y' [coordinates of the point snapped to
 the line], 'd' [distance from the input point to the new snapped
 point], and 'chain' [accumulated along-line distance from the starting
 of the line to the snapped point]
}
\author{Javier Garcia-Pintado}
\seealso{\code{\link[sp]{Spatial-class}}}
\keyword{ utilities }
\keyword{ spatial }
