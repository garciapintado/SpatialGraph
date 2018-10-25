\name{sgChVIDs}
\alias{sgChVIDs}
\title{Change vertex IDs in a SpatialGraph}
\description{Change the field "ID" in the vertex slot, \code{v}, of a
  \link{SpatialGraph}. The fields \code{v0} and \code{v1} of the edge
  slot, \code{e}, are accordingly updated
}
\usage{
  sgChVIDs(obj, IDa, IDp = NULL)
}
\arguments{
  \item{obj}{A \link{SpatialGraph} object}
  \item{IDa}{A vector indicating the updated vertex IDs}
  \item{IDp}{A vector indicating the prior vertex IDs}
}
\details{
  If \code{IDp} is not provided, it is assumed that the vector of
  updated indexes is sorted equally to the order in which the vertices
  are stored in the slot \code{v} of the \link{SpatialGraph}. If
  \code{IDp} is provided, the mapping IDp -> IDa is used
  for reclassifying the vertices.
}
\value{
  A \link{SpatialGraph} object
}
\keyword{ utilities }
\keyword{ spatial }
