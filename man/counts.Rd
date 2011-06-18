\name{counts}
\docType{methods}
\alias{counts,ExonCountSet-method}
\alias{counts<-,ExonCountSet,matrix-method}
\title{
   Accessors for the 'counts' slot of a ExonCountSet object.
}
\description{
   The counts slot holds the count data as a matrix of non-negative integer
   count values, one row for each observational unit (an exon of a gene), and
   one column for each sample.
}
\usage{
\S4method{counts}{ExonCountSet}(cds, normalized=FALSE)
\S4method{counts}{ExonCountSet,matrix}(cds) <- value
}
\arguments{
  \item{cds}{an \code{ExonCountSet} object.}
  \item{normalized}{logical indicating whether or not to divide the
    counts by the size factors before returning.}
  \item{value}{an integer matrix of counts, each row corresponding to an exon and each column corresponding to a sample.}
}
\examples{
   data("pasillaExons", package="pasilla")
   head( counts( pasillaExons ) )
}   
