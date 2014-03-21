\docType{class}
\name{DEXSeqResults-class}
\alias{DEXSeqResults}
\alias{DEXSeqResults-class}
\title{DEXSeqResults object}
\usage{
  DEXSeqResults( object )
}
\arguments{
  \item{object}{
    A \code{DEXSeqDataSet} object.
  }
}
\value{
  A DEXSeqResults object.
}
\description{
  The \code{DEXSeqResults} object is a subclass 
  of a \code{DataFrame}. It collects relevant 
  information from a \code{DEXSeqDataSet} object with
  the results generated from testing for differences
  in exon usage.
}
\examples{
  data(pasillaDEXSeqDataSet, package="pasilla")
  dxd <- estimateSizeFactors( dxd )
  dxd <- estimateDispersions( dxd )
  dxd <- testForDEU( dxd )
  dxr <- DEXSeqResults( dxd )
}
