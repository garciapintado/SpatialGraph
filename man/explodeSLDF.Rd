\name{explodeSLDF}
\alias{explodeSLDF}
\title{Explode Lines in a SpatialLinesDataFrame}
\description{explode Lines in a SpatialLinesDataFrame, so that each
  single Line, within each Lines slot, is upgraded as a new 1-Line Lines slot
}
\usage{
explodeSLDF(SLDF, FID)
}
\arguments{
  \item{SLDF}{a SpatialLinesDataFrame}
  \item{FID}{if not NULL, field name, within the attribute table considered as additional unique identifier, so that incremental numeric values will added to this field to avoid duplicate values}
}
\value{
  a SpatialLinesDataFrame
}
\keyword{ utilities }
\keyword{ spatial }
