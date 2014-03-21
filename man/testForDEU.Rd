\name{testForDEU}
\alias{testForDEU}
\title{Test for Differential Exon Usage.}
\description{
  This will perform a likelihood ratio test for differential 
  exon usage. Internally, it calls the DESeq2 function 
  \code{nbinomLRT}. 
}
\usage{

  testForDEU( object, 
     fullModel = design(object), 
     reducedModel = ~ sample + exon, 
     BPPARAM=MulticoreParam(workers=1) )

}
\arguments{
  \item{object}{
     A DEXSeqDataSet object.
  }
  \item{fullModel}{
     The full model formula
  }
  \item{reducedModel}{
     Null model formula.
  }
  \item{BPPARAM}{
     A "BiocParallelParam" instance. 
     See \code{?bplapply} for details.
  }
}
\details{

   The information of the variables of the formulas 
   should be present in the \code{colData} of the 
   \code{DEXSeqDataSet} object.

}
\value{

   A \code{DEXSeqDataSet} with slots filled with information
   about the test.

}
\examples{

  data(pasillaDEXSeqDataSet, package="pasilla")
  dxd <- estimateSizeFactors( dxd )
  dxd <- estimateDispersions( dxd )
  dxd <- testForDEU( dxd )

}

