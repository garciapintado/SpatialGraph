\name{attSGe}
\alias{attSGe}
\title{Add or Modify attributes in SpatialGraph edges}
\description{Add or Modify attributes in SpatialGraph edges
}
\usage{
attSGe(SG, att, eID, val, default)
}
\arguments{
  \item{SG}{\link{SpatialGraph}}
  \item{att}{name of the field [column] in the edge dataframe to be added/modified}
  \item{eID}{edge identifiers [row.names of the edge data.frame]}
  \item{val}{values corresponding the \code{eID} above}
  \item{default}{default values for edges not considered in \code{eID} above}
}
\value{
  A \link{SpatialGraph}
}
\keyword{ utilities }
\keyword{ spatial }
