\name{distSG}
\alias{distSG}
\title{Calculate across-network distance for a set of sparse points}
\description{This function obtains the across-network distance for a set
  of sparse points, by using the distance slot in a \link{SpatialGraph}. The
  calculation is supported by a previously calculated between vertex
  distance matrix
  [via a call to the library \code{igraph} by the function \link{distSGv}].
  The SpatialGraph is considered as undirected for distance
  calculation.
  If \code{euc=TRUE} [default], the distance between two points is defined within this function
  as the maximum of both the minimum along-network distance and the Euclidean
  distance. The distance itself between the points
  in x,y and the network is neglected in the function for the
  along-network distance.
  Both, \code{x} and \code{y}, are \code{SpatialPointsDataFrame}
  objects, which must contain at least the fields
  \code{eID} and \code{chain}, which describe their relationship with
  the \code{SpatialGraph} object defined by \code{SG}. These can be
  obtained with either the function \code{\link{pointsSLDFchain}} or
  \link{pointsToLines} (the latter is faster, but depends on GEOS)
}
\usage{
 distSG(SG, x, y = NULL, euc = TRUE, wei = NULL, getpath = FALSE)
}
\arguments{
  \item{SG}{\link{SpatialGraph}}
  \item{x}{\code{SpatialPointsDataFrame}}
  \item{y}{\code{SpatialPointsDataFrame}}
  \item{euc}{boolean scalar, whether to use Euclidean distance as
    minimum threshold for resulting distances}
  \item{wei}{if not null, field in \code{SG@e} with a variable to obtain
    a state-related weight. See details below.}
  \item{getpath}{if TRUE (and wei != NULL), \code{eID} identifiers for
    each path from \code{x} to \code{y} elements is returned}
}
\value{
  If \code{wei=NULL}, a matrix of distances between \code{x} and
  \code{y}. If \code{wei} is not \code{NULL}, a list with a distance
  matrix and weight matrix (plus a matrix with {eID} identifiers
  for the path, if \code{getpath=TRUE}) is returned.
}
\details{
The application of state-related weights in this version is a simple
state-dependent weight matrix related to some field in \code{SG@e} [i.e. the
edges in the input \code{SpatialGraph}]. The only current calculation
evaluates the path between queried points (x,y), and along the path, for
every junction and jump into a new edge, the ratio for the evaluated state variable (taken as
the highest value divided by the lowest value) between the two edges at
the junction is obtained. Currently a
maximum ratio equal to 10.0 is hard-coded. The product
of ratios along the path gives the weight.
}
\examples{
  if (1 > 2) { # not run
    dem <- readGDAL(file.path(system.file('external',package='hydrosim'),
                    'watershed1','IDRISI_maps','dem','dem.rst'))   # SpatialGridDataFrame
    plotGmeta(layer=dem, xlim=662500 + 2500 * c(-1,+1),
              ylim=4227500 + 2500 * c(-1,1), zlim='strloc', as.na=0)

    # generate some crossing lines
    zz <- list()
    zz[[1]] <- digitGmeta(layer=dem, type='Lines', ID=1)
    zz[[2]] <- digitGmeta(layer=dem, type='Lines', ID=2)
    zz[[3]] <- digitGmeta(layer=dem, type='Lines', ID=3)
    SL <- SpatialLines(zz)
    SG <- sl2sg(SL, getpath=TRUE)
    points(SG@v, cex=2)                    # plot SpatialGraph vertices

    apath <- SG@path[[1,2]]                # iteratively plot a path as an example
    for (iv in 1:length(apath$v)) {
      points(SG@v[apath$v[iv],], cex=2,pch=2)
      if (iv == length(apath$v))
        break
      lines(SG@e[apath$e[iv],],col='blue',lwd=2,lty=2)
      Sys.sleep(1)
    }

    # sample a few points [as a matrix] close to some edges
    xy    <- digit()                    # sample locations
    xych  <- pointsToLines(xy, SG@e)    # SpatialPointsDataFrame mapping
    points(xy, col='blue', pch=3)
    points(xych, col='darkgreen', pch=19)

    # along-network distance
    xyndis <- distSG(SG, xych)

    # state-dependent weighted along-network distance
    SG@e@data$wxs <- 3+round(runif(nrow(SG@e@data)),2)         # [m2] foo wetted cross-section areas
    SG@e@data

    xywdis <- distSG(SG, xych, wei='wxs')
    xywdis <- xywdis$dis * xywdis$wei       # Schur weight application into distance estimation
  }
}
\author{Javier Garcia-Pintado, e-mail: \email{jgarciapintado@marum.de}}
\keyword{ utilities }
\keyword{ spatial }
