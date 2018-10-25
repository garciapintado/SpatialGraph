\name{SpatialGraph-class}
\docType{class}
\alias{SpatialGraph-class}

\title{Class "SpatialGraph"}
\description{
  Class for spatial networks
}
\section{Objects from the Class}{
Objects can be created by calls to the function \link{SpatialGraph}
}
\section{Slots}{
  \describe{
    \item{\code{v}:}{Object of class \code{"SpatialPointsDataFrame"},
      whose data.frame must contain the "ID" field as unique identifier}
    \item{\code{e}:}{Object of class \code{"SpatialLinesDataFrame"},
      whose data.frame must contain the fields \code{v0} and \code{v1}
      matching the unique identifiers "ID" in the slot \code{v} data.frame}
    \item{\code{dist}:}{Matrix, representing the undirected along-graph distance
      between all vertices in the network}
    \item{\code{path}:}{list with variable length arrays describing the
      minimum distance path between vertices}


  }
}
\author{Javier Garcia-Pintado, e-mail:
  \email{j.garcia-pintado@reading.ac.uk}}
\keyword{ classes }
\keyword{ spatial }
