\name{pointLineD}
\alias{pointLineD}
\title{Euclidean distance from a set of points to a line segment}
\description{
  \code{pointLineD} returns a list with a number of components from a points to line segment analysis}
\usage{
pointLineD(xy, xyp)
}
\arguments{
  \item{xy}{2 x 2 [x,y] matrix defining the start and end of the segment}
  \item{xyp}{p x 2 [x,y] matrix with a point set}
}
\details{
  \code{pointLineD} conduct a detailed points to segment distance analysis, returned as a list
}
\value{
 A list with the input components \code{xy} and \code{xyp}, and the aditional components: \code{d}, point-line distance (distance between the points in \code{xyp} and their perpendicular projections of the line); \code{dc},  diferential chainage over [x0,y0] (> 0 if the projection goes in the segment direction); \code{cross}, boolean vector indicating whether the perpendicular projection of the points crosses the segment, or not
}
\author{Javier Garcia-Pintado}
\seealso{\code{\link[sp]{Spatial-class}}}
\keyword{ utilities }
\keyword{ spatial }
