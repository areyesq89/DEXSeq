\name{perGeneQValue}
\alias{perGeneQValue}
\title{Summarize per-exon p-values into per-gene q-values.}
\description{
  The use case for this function is the following analysis: given
  per-exon p-values for null hypothesis H0, we can determine the number
  of genes in which at least for one exon H0 is rejected. What is the
  associated false disovery rate?
}
\usage{
perGeneQValue(object, p = "pvalue", method = perGeneQValueExact)
}
\arguments{
  \item{object}{
    An \code{DEXSeqResults} object.
  }
  \item{p}{A character string indicating the name of the columns in
    \code{DEXSeqResults} from which to take the per-exon p-values.}
  \item{method}{Use the default value. This is for debugging only.}
  }
\details{Details}
\value{
  A named numeric vector, values are per-gene q-values, names are gene.
}
\examples{

  data(pasillaDEXSeqDataSet, package="pasilla")
  dxd <- estimateSizeFactors( dxd )
  dxd <- estimateDispersions( dxd )
  dxd <- testForDEU( dxd )
  dxr <- DEXSeqResults( dxd )

  perGeneQValue(dxr)


}
