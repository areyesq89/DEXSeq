\name{plotMA-methods}
\docType{methods}
\alias{plotMA}
\alias{plotMA-methods}
\alias{plotMA,DEXSeqDataSet-method}
\alias{plotMA,DEXSeqResults-method}
\title{Generate an MA plot}
\description{
  This function generates an MA plot.
}
\usage{

  \S4method{plotMA}{DEXSeqDataSet}( object, alpha=0.1, ylim=c(-2, 2), foldChangeColumn=NULL, ... )

  \S4method{plotMA}{DEXSeqResults}( object, alpha=0.1, ylim=c(-2, 2), foldChangeColumn=NULL, ... )

}
\arguments{
  \item{object}{
     Either a \code{DEXSeqDataSeq} or a \code{DEXSeqResults}.
  }
  \item{alpha}{
     the false discovery rate, i.e., threshold to the adjusted p values, to be used 
     to colour the dots
  }
  \item{ylim}{
     y-limits of the plot
  }
  \item{foldChangeColumn}{
     Name of the column containing the fold changes to plot. If NULL,
     the first fold change column from the object will be taken.
  }
  \item{...}{
    Further parameters to be passed through to \code{\link{plot}}.
  }
}

\examples{
 \dontrun{
    data(pasillaDEXSeqDataSet, package="pasilla")
    dxd <- estimateSizeFactors( dxd )
    dxd <- estimateDispersions( dxd )
    dxd <- testForDEU( dxd )
    dxd <- estimateExonFoldChanges( dxd )
    dxr <- DEXSeqResults( dxd )

    plotMA( dxr )

    plotMA( dxd )
    
 }
}

