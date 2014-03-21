\docType{methods}
\name{counts}
\alias{counts}
\alias{counts,DEXSeqResults-method}
\alias{counts,ExonCountSet-method}
\alias{counts,ExonCountSet-method,deprecated}
\title{Accessors for the 'counts' slot of a DEXSeqResults object.}
\usage{
\S4method{counts}{DEXSeqResults}(object,normalized=FALSE)
}
\arguments{
  \item{object}{a \code{DEXSeqResults} object.}

  \item{normalized}{logical indicating whether or not to
  divide the counts by the size factors or normalization
  factors before returning (normalization factors always
  preempt size factors)}
}
\value{
  an integer matrix
}
\description{
The counts slot holds the count data as a matrix of
non-negative integer count values, one row for each
observational unit (gene or the like), and one column for
each sample.
}
\examples{
  data(pasillaDEXSeqDataSet, package="pasilla")
  head( counts( dxd ))
}
