\name{pointPolylineD}
\alias{pointPolylineD}
\title{closest points in a polyline to a set of points}
\description{
  \code{pointPolylineD} returns a list with a number of components from a points to polyline analysis}
\usage{
pointPolylineD(xy, xyp)
}
\arguments{
  \item{xy}{n x 2 [x,y] matrix defining the polyline}
  \item{xyp}{p x 2 [x,y] matrix with a point set}
}
\details{
  \code{pointPolylineD} conducts a detailed points to polyline distance analysis. First the distance from the set of points to the lines defined by every single segment in the polyline is obtained by succesive calls to \code{pointLineD}, then the distance to every single node in the polyline are also obtained. The lower distance is chosen.
}
\value{
 A \code{data.frame} with the columns: \code{inode} is the index of the first node in the closest segment to each point, \code{x0} and \code{y0} are the corresponding coordinates of those nodes, \code{xc} and \code{yc} are the coordinates of the point in the polyline closest to each point in \code{xyp}, these may be but are not necessarily one the polyline nodes, \code{dis} it the distance from each point tho the polyline, \code{chain0} is the chainage of \code{x0,y0} with the polyline, and \code{dc} is the differential chainage from \code{xc,yc} to \code{x0,y0}
}
\seealso{\code{\link[sp]{Spatial-class}}}
\keyword{ utilities }
\keyword{ spatial }
