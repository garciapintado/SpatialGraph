\name{polylineChainage}
\alias{polylineChainage}
\title{Obtain the chainage of nodes along a polyline}
\description{
  Obtain the chainage of nodes along a polyline [2-col matrix]
}
\usage{
polylineChainage(xy)
}
\arguments{
  \item{xy}{a 2-column matrix representing the polyline nodes}
}
\details{
 \code{polylineChainage} calculates a vector of chainage values
 [along-polyline distances] from each node in a polyline to the initial node 
}
\value{
  A vector
}
\author{Javier Garcia-Pintado, e-mail:
  \email{j.garcia-pintado@reading.ac.uk}}
\seealso{\code{\link{polylineLength}}}
\keyword{ utilities }
\keyword{ spatial }
