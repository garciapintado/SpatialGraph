\name{distSGv}
\alias{distSGv}
\title{Calculate the distance slot in a SpatialGraph}
\description{Calculate the distance slot in a \link{SpatialGraph}. This
  is done via a call to the library \code{igraph}, which does the
  calculation. Distances are undirected.}
\usage{
 distSGv(SG, getpath = FALSE)
}
\arguments{
  \item{SG}{\link{SpatialGraph}}
  \item{getpath}{boolean. Whether to calculate the \code{SG@path} slot}
}
\value{
  A \link{SpatialGraph} with the slot \code{dist} (and
  \code{path} if requested) recalculated
}
\keyword{ utilities }
\keyword{ spatial }
