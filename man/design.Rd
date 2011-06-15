\name{design}
\alias{design}
\alias{design<-}
\title{
    Accessor function for the design information in a ExonCountSet object
}
\description{
    The design vector is a factor that assigns to each column of the count data
    a condition (or treatment, or phenotype, or the like). This information is stored
    in the ExonCountSet's "phenoData" slot as a row named "design".
}
\usage{
design(ecs, drop=TRUE, asAnnotatedDataFrame=FALSE)
design(ecs) <- value
}
\arguments{
  \item{ecs}{
      An ExonCountSet
   }
  \item{drop}{
      This option gives the option that if it is a single factor design, the value can 
      be presented as a data frame instead of a single factor value.
   }
  \item{asAnnotatedDataFrame}{
      Indicates whether the result should be presented as an AnnotatedDataFrame.
   }
  \item{value}{
      A vector or matrix with conditions for the samples, one row for each column in the count data.
   }
}
\author{
   Simon Anders, sanders@fs.tum.de
}
\examples{
   library(DEXSeq)
   data("pasillaExons", package="pasilla")
   design( pasillaExons )
}
