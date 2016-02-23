\docType{class}
\name{DEXSeqResults-class}
\alias{DEXSeqResults}
\alias{DEXSeqResults-class}
\title{DEXSeqResults object}
\usage{
  DEXSeqResults( object, independentFiltering=TRUE, filter)
}
\arguments{
    \item{object}{
        A \code{DEXSeqDataSet} object.
    }
    \item{independentFiltering}{
        Logical indicating whether independent filtering should be 
    	applied automatically. 
	}
    \item{filter}{
        A vector of filter statistics over which the independent
        filtering will be optimized. The default is the normalized
        exon means.
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
