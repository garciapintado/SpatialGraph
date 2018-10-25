\name{splitSLDF}
\alias{splitSLDF}
\title{Split 1-Line Lines in a SpatialLines or a SpatialLinesDataFrame by intersection
  with a point dataset}
\description{
  \code{splitSLDF} divides the 1-Line Lines in the \code{SpatialLines} or the
  \code{SpatialLinesDataFrame} at intersections with the input point dataset
}
\usage{
splitSLDF(SLDF, SPDF, dmax=NULL)
}
\arguments{
  \item{SLDF}{length-1 SpatialLinesDataFrame or SpatialLines object}
  \item{SPDF}{SpatialPointsDataFrame}
  \item{dmax}{maximum distance between points in \code{SPDF} and the
    polylines in \code{SLDF}, for these to be considered for poyline splitting}
}
\details{
  \code{splitPolyline} obtain the closest points in the SpatialLinesDataFrame to a
  given input set of points. Those closest points are used to divide the
  polylines in a number of transects. The individual transects are
  clipped to the input point dataset, so the different transects are
  continuous in space. Note that if the input points is quite appart
  from the polyline, the output sequence of transects may substantially
  differ form the input polyline at rupture zones. The input parameter
  \code{dmax} is provided as a mean to avoid too strange splitting
  results. Setting dmax to a ver low value will reduce the spureous
  results, but also the input points need to be closer to the lines for
  the adequate recognition of splitting points
}
\value{
 A SpatialLinesDataFrame or a SpatialLines, according to the input
}
\author{Javier Garcia-Pintado}
\seealso{\code{\link[sp]{Spatial-class}}}
\keyword{ utilities }
\keyword{ spatial }
