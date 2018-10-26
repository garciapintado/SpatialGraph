\name{textSGe}
\alias{textSGe}
\title{Label edges in a SpatialGraph plot}
\description{
  A \link{SpatialGraph} contains a \code{SpatialLinesDataFrame},
  describing the network topology. This function adds line \code{IDs}
  and direction arrows to an existing plot of a \link{SpatialGraph}.
}
\usage{
textSGe(SG, acol='wheat', tcol='navyblue', arr.length=0.4)
}
\arguments{
  \item{SG}{\code{SpatialGraph}}
  \item{acol}{color of the graph direction arrows}  
  \item{tcol}{color of the text for graph edge IDs}
  \item{arr.length}{length of the direction arrows}
}
\value{
 Arrows and edge \code{IDs} added to a \link{SpatialGraph} plot
}
\keyword{ utilities }
\keyword{ spatial }
