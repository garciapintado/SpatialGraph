\name{pointsToLines}
\alias{pointsToLines}
\title{Snap a set of points to a set of lines}
\description{
  This function snaps a set of points to a set of lines based on the
  minimum distance of each point to any of the lines}
\usage{
pointsToLines(points, lines, withAttrs = TRUE, withDis = TRUE, withChain = TRUE)
}
\arguments{
  \item{points}{An object of the class SpatialPoints or
    SpatialPointsDataFrame, or a 2-col matrix of [x,y] coordinates}
  \item{lines}{An object of the class SpatialLines or
    SpatialLinesDataFrame}
  \item{withAttrs}{Boolean value for preserving (TRUE) or getting rid
    (FALSE) of the original point attributes. Default: TRUE. This
    parameter is optional}
  \item{withDis}{Boolean value for including distance from source points
  to snapped-to-lines points}
  \item{withChain}{Boolean value for including the chainage of the
    snapped points in their corresponding lines}
}
\value{
 A SpatialPointsDataFrame object as defined by the R package 'sp'. This
 object contains the snapped points, therefore all of them lie on the
 lines. The returned object contains the fields 'lid', 'eID', and
 'chain', providing information about the relationship between the
 source data points, the snapped data points, and its location within
 the network: 'lid', and 'eID' are the line index and line ID,
 respectively, of the lines in which the new snapped points lie; 'dis'
 is the distance between the input points and the snapped points, and 'chain'
 is the chainage of the snapped point within the corresponding line
}
\author{Javier Garcia-Pintado}
\seealso{\code{\link[sp]{Spatial-class}}}
\keyword{ utilities }
\keyword{ spatial }
