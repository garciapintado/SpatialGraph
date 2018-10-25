\name{routeSDG}
\alias{routeSDG}
\title{Accumulate sources/sinks along a directed SpatialGraph}
\description{
 Assume a SpatialGraph is directed and conduct an accumulation of
 source/sink values at nodes across the network. The accumulation
 assumes no delay in transmission 
}
\usage{
routeSDG(SDG, FUN='cumsum', ifld='inflow')
}
\arguments{
  \item{SDG}{\code{SpatialGraph}, assumed as directed}
  \item{FUN}{name of a function to be applied for the routing}
  \item{ifld}{name on the field in the \code{SpatialPointDataFrame} vertex slot to be used used as source/sink}
}
\details{
  The \code{SpatialGraph}, used as input, must have the \code{ifld} field
  to be used as input, in the vertices slot \code{v} (a SpatialPointsDataFrame). The accumulated output is provided
  as the new field \code{ofld} in \code{v}. The edges slot \code{e}
  serves to route the input across the network
}
\value{
A \code{SpatialGraph} with the added \code{ofld} field in the vertex slot
}
\keyword{ utilities }
\keyword{ spatial }
