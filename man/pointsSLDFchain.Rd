\name{pointsSLDFchain}
\alias{pointsSLDFchain}
\title{Obtain chainage from sparse points along a SpatialLinesDataFrame}
\description{For a set of points, obtains the closest Line object in a
  SpatialLinesDataFrame. The function assumes that each Feature (entry
  in the DataFrame part of the SpatialLinesDataFrame) just contains one
  \code{Line} (i.e. one polyline). The within-polyine chainage (that is,
  distance from the initial point of the poyline to the mapping of the
  point into the polyline) is also returned. If mask is NULL, each point
  in the set is assigned a line in SLDF by Euclidean distance. If mask
  is provided, the match between \code{mask} and the \code{SLmsk}
  field in \code{SLDF} is used instead for polyline assignation.
}
\usage{
 pointsSLDFchain(SLDF, xy, SLmsk='FEAT_ID', mask=NULL,  type='SpatialPointsDataFrame')
}
\arguments{
  \item{SLDF}{\code{SpatialLinesDataFrame}}
  \item{xy}{REAL [n,2] matrix of points, or a \code{SpatialPointsDataFrame}}
  \item{SLmsk}{is !is.null(mask) this is the field in the SLDF data.frame matching the values in mask}
  \item{mask}{REAL, OPT, [n] a vector indicating to which line in SLDF is
	      related each point}
  \item{type}{character. Either 'SpatialPointsDataFrame' or
    'mapping'. In the latter case, just the chainage in line feature
    identifiers are returned}
}
\value{
  A data.frame with two columns, 'chai', and 'eIDs', where 'eIDs' are
  the row names of the data.frame component of the input SpatialLinesDataFrame
}
\author{Javier Garcia-Pintado, e-mail: \email{jgarciapintado@marum.de}}
\keyword{ utilities }
\keyword{ spatial }
