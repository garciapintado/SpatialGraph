\name{sg2igraph}
\alias{sg2igraph}
\title{Map a SpatialGraph into an igraph}
\description{
  The vertex and edge information in  a \code{SpatialGraph} is
  mapped into an \code{igraph} object
}
\usage{
sg2igraph(sg, directed=FALSE)
}
\arguments{
  \item{sg}{\code{SpatialGraph}}
  \item{directed}{whether the resulting \code{igraph} is directed}
}
\details{
  It is assumed that the \code{SpatialGraph}, used as input, is correct
  (i.e.g all records in \code{sg@e@data} have the two first field correctly
  identifying the field 'ID' in \code{sg@v}. It is also assumed that the
  \code{sg@e@data} data.frame has the fields \code{div} and
  \code{len}. These two are highly useful to conduct network operations
  on the resulting \code{igraph}
}
\value{
An \code{igraph}
}
\keyword{ utilities }
\keyword{ spatial }
