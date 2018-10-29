\name{sl2sg}
\alias{sl2sg}
\title{map a SpatialLinesDataFrame into a SpatialGraph}
\description{
  This function is the major workhorse to map an input
  \code{SpatialLinesDataFrame}, as defined in the package sp, into a
  \code{SpatialGraph} by using the spatial connectivity. Input is first
  exploded by using \code{\link{explodeSLDF}}, and then all vertices in the
  \code{SpatialGraph} are automatically generated according to crossings
  in the input polylines.
}
\usage{
sl2sg(SL, clipd = NULL, getdist = TRUE, getpath = FALSE)
}
\arguments{
  \item{SL}{\code{SpatialLinesDataFrame} as defined in package sp}
  \item{clipd}{distance threshold for clipping features, If NULL, a
    value of 1.0E-04 of the domain side size is used}
  \item{getdist}{calculate the \code{dist} slot in the returned
    \code{SpatialGraph}}
 \item{getpath}{calculate the \code{path} slot in the returned
    \code{SpatialGraph}}
}
\details{
  A \code{SpatialGraph} is generated
}
\value{
A \code{SpatialGraph}
}
\examples{
  #  x   y
  # create list of Line objects
 if (1 > 2) {
  library(sp)
  library(SpatialGraph)
  zz <- list()
  zz[[1]] <- Line(matrix(
   c(661750, 4229150,
     662650, 4229450,
     663550, 4227650,
     663550, 4226850), ncol=2, byrow=TRUE))
  zz[[2]] <- Line(matrix(
   c(660250, 4229650,
     661050, 4226450,
     662550, 4225350,
     664850, 4225850,
     664650, 4229150,
     662350, 4228850), ncol=2, byrow=TRUE))
  # upgrade Line as Lines
  for (i in 1:length(zz)) {
    zz[[i]] <- Lines(list(zz[[i]]), ID=i)
  }
  # as SpatialLines
  SL <- sp::SpatialLines(zz)
  # as SpatialGraph including path calculation
  SG <- sl2sg(SL, getpath=TRUE) 

  plot(SL, axes=TRUE)  
  points(SG@v, cex=2)
  lines(SG@e, lwd=2)
  points(SG@v, cex=2, col='grey', pch=19)
  text(SG@v, labels=SG@v$ID)
  # label edges and directions
  textSGe(SG)
  # show a distance matrix between nodes
  SG@dist
  # show path from node 1 to 3
  SG@path[1,3]
 }
}
\author{Javier Garcia-Pintado, e-mail:
  \email{j.garcia-pintado@marum.de}}
\keyword{ utilities }
\keyword{ spatial }
