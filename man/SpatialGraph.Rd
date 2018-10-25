\name{SpatialGraph}
\alias{SpatialGraph}
\title{Create a SpatialGraph object}
\description{
  A \code{SpatialGraph} object is created
}
\usage{
SpatialGraph(v, e, dist = NULL, path = NULL)
}
\arguments{
  \item{v}{\code{SpatialPointsDataFrame}}
  \item{e}{\code{SpatialLinesDataFrame}}
  \item{dist}{along-network (symmetric) distance matrix}
  \item{path}{matrix of lists with paths corresponding to
    \code{dist}. While distances between
  vertex couples are symmetric, the path matrix is not symmetric as
  individual path to from source vertex to destination vertex. Each list in the matrix has two S3
  components (v,e) describing vertices (including bounds) and edges
  along the path. Thus it is always one less edge than then number of
  vertices in the path}
}
\value{
  \code{SpatialGraph} returns an object of class \link{SpatialGraph-class}
}
\keyword{ classes }
\keyword{ spatial }
