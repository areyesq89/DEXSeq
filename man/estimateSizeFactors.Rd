\docType{methods}
\name{estimateSizeFactors}
\alias{estimateSizeFactors}
\alias{estimateSizeFactors,DEXSeqDataSet-method}
\alias{estimateSizeFactors,ExonCountSet-method}
\alias{estimateSizeFactors,ExonCountSet-method,deprecated}
\title{Estimate the size factors for a DEXSeqDataSet}
\usage{
  \S4method{estimateSizeFactors}{DEXSeqDataSet}(object,locfunc=median,geoMeans)
}
\arguments{
  \item{object}{a DEXSeqDataSet}

  \item{locfunc}{a function to compute a location for a
  sample. By default, the median is used. However,
  especially for low counts, the
  \code{\link[genefilter]{shorth}} function from the
  genefilter package may give better results.}

  \item{geoMeans}{by default this is not provided and the
  geometric means of the counts are calculated within the
  function. A vector of geometric means from another count
  matrix can be provided for a "frozen" size factor
  calculation}
}
\value{
  The DEXSeqDataSet passed as parameters, with the size
  factors filled in.
}
\description{
  Estimate the size factors for a DEXSeqDataSet
}
\details{
  This function calls the method estimateSizeFactors
  for the DESeqDataSet object, and adapts it to the
  specific use to the DEXSeqDataSet.
}
\references{
  \itemize{ \item Simon Anders, Wolfgang Huber: Differential
  expression analysis for sequence count data. Genome Biology
  11 (2010) R106,
  \url{http://dx.doi.org/10.1186/gb-2010-11-10-r106} }
}

