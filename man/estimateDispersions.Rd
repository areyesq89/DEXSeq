\docType{methods}
\name{estimateDispersions}
\alias{estimateDispersions}
\alias{estimateDispersions,DEXSeqDataSet-method}
\title{Estimate the dispersions for a DEXSeqDataSet}
\usage{

\S4method{estimateDispersions}{DEXSeqDataSet}( object, fitType=c("parametric","local","mean"), maxit=100, niter=10, quiet=FALSE, formula=design(object), BPPARAM=SerialParam())

}
\arguments{
  \item{object}{
    A DEXSeqDataSet
  }
  \item{fitType}{
    Either "parametric", "local", or "mean"
    for the type of fitting of dispersions to the mean
    intensity. See ?estimateDispersions,DESeqDataSet-method
    for details.
  }
  \item{maxit}{
    Control parameter: maximum number of 
    iterations to allow for convergence
  }
  \item{niter}{
    Number of times to iterate between estimation of means
    and estimation of dispersion.
  }
  \item{quiet}{
    Whether to print messages at each step
  }
  \item{formula}{ 
    Formula used to fit the dispersion estimates 
  }
  \item{BPPARAM}{ 
     A "BiocParallelParam" instance. 
     See \code{?bplapply} for details.
  }
}
\value{
  A DEXSeqDataSet with the dispersion information 
  filled in as metadata columns.
}
\description{
  This function obtains dispersion estimates for negative
  binomial distributed data for the specific case for DEXSeq.
}
\details{
  See ?estimateDispersions,DESeqDataSet-method for details.
}
\examples{
  data(pasillaDEXSeqDataSet, package="pasilla")
  dxd <- estimateSizeFactors( dxd )
  dxd <- estimateDispersions( dxd )
}

