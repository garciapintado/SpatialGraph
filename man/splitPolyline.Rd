\name{splitPolyline}
\alias{splitPolyline}
\title{Split a polyline into a number of transects}
\description{
  \code{splitPolyline} returns a list with a number of transects along a
polyline}
\usage{
splitPolyline(xy, xyp, dmax)
}
\arguments{
  \item{xy}{2-column [x,y] matrix defining the polyline nodes}
  \item{xyp}{2-column [x,y] matrix with a point set}
  \item{dmax}{maximum distance between points in \code{xy} and the
    polyline, for these to be considered for poyline splitting}
}
\details{
  \code{splitPolyline} obtain the closest points in a polyline to a
  given input set of points. Those closest points are used to divide the
  polyline in a number of transects. The indivudual transects are
  clipped to the input point dataset, so the different transects are
  continuous in space. Note that if the input points is quite appart
  from the polyline, the output seqence of transect may substantially
  differ form the input polyline at rupture zones
}
\value{
 A list in which each element is a matrix representing an individual polyline
}
\author{Javier Garcia-Pintado}
\seealso{\code{\link[sp]{Spatial-class}}}
\keyword{ utilities }
\keyword{ spatial }
