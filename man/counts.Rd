\name{counts}
\alias{counts}
\alias{counts<-}
\title{
   Accessors for the 'counts' slot of a ExonCountSet object.
}
\description{
   The counts slot holds the count data as a matrix of non-negative integer
   count values, one row for each observational unit (an exon of a gene), and
   one column for each sample.
}
\usage{
counts(ecs, normalized=FALSE)
counts(ecs) <- value
}
\arguments{
  \item{ecs}{
     An ExonCountSet object
   }
  \item{normalized}{
     Logical indication if the counts want to be normalized by the size factor.
   }
  \item{value}{
     A matrix containing integer counts, each row corresponding to an exon and each column corresponding to a sample.
   }
}
\examples{
   library(DEXSeq)
   data("pasillaExons", package="pasilla")
   head( counts( pasillaExons ) )
}   
