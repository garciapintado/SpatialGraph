\name{revSGe}
\alias{revSGe}
\title{Reverse Lines in a SpatialGraph}
\description{
  A \link{SpatialGraph} contains a \code{SpatialLinesDataFrame},
  describing the network topology. The input \code{eID} indicates the
  identifiers of a set of lines (edges) in the network to be reversed. Note
  \code{eID} does not refer to the line index within \code{SG@e}, but to
  the Feature Identifiers, as extracted from \code{row.names(SG@e@data)}
}
\usage{
revSGe(SG, eID)
}
\arguments{
  \item{SG}{\code{SpatialGraph}}
  \item{eID}{vector of Feature Identifiers for lines to be reversed}
}
\details{
  Note \code{eID} does not refer to the line index within \code{SG@e}, but to
  the Feature Identifiers, as extracted from
  \code{row.names(SG@e@data)}. Accordingly to the reversed coordiantes,
  the corresponding fields ["v0","v1"], are interchanged.
}
\value{
 A \link{SpatialGraph}
}
\keyword{ utilities }
\keyword{ spatial }
