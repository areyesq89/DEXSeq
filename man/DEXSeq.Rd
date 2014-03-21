\name{DEXSeq}
\alias{DEXSeq}
\title{Performs the differential exon usage test in a single command}
\description{
  This function is a wrapper that calls the necessary functions
  to perform a differential exon usage test in a single command.
}
\usage{
  DEXSeq(object, 
    fullModel=design(object),
    reducedModel = ~ sample + exon, 
    BPPARAM=MulticoreParam(workers=1), 
    fitExpToVar="condition", quiet=TRUE )
}
\arguments{
  \item{object}{
    An DEXSeqDataSet object.
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
  \item{fitExpToVar}{
    A variable name contained in the sample data. The expression 
    values will be fitted to this variable using the the formula 
    " ~ sample + fitExpToVar * exon".
  }
  \item{quiet}{
    Whether to print messages at each step
  }
}
\value{
  A \code{DEXSeqResults} object.
}
\examples{
  data(pasillaDEXSeqDataSet, package="pasilla")
  dxr <- DEXSeq( dxd )
}
