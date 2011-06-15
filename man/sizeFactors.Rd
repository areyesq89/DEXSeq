\name{sizeFactors}
\alias{sizeFactors}
\alias{sizeFactors<-}
\title{
    Accessor functions for the sizeFactors information in a ExonCountSet
}
\description{
    The sizeFactors vector assigns to each column of the count data
    a value, the size factor, such that count values in the columns can be brought to
    a common scale by dividing by the corresponding size factor.
}
\usage{
sizeFactors(ecs)
sizeFactors(ecs) <- value
}
\arguments{
  \item{ecs}{
      An ExonCountSet
}
  \item{value}{
      a vector of number, one size factor for each column in the count data
}
}

\seealso{
\code{\link{estimateSizeFactors}}
}
\examples{
   data("pasillaExons", package="pasilla")
   pasillaExons <- estimateSizeFactors( pasillaExons )
   sizeFactors(pasillaExons)
}
\author{
   Simon Anders, sanders@fs.tum.de
}

