\name{plotMA-methods}
\docType{methods}
\alias{plotMA}
\alias{plotMA-methods}
\alias{plotMA,data.frame-method}
\alias{plotMA,ExonCountSet-method}
\title{Generate an MA plot}

\description{
      This function generates an MA plot.
}
\usage{

\S4method{plotMA}{data.frame}( object, ylim = NULL,
  colNonSig = "gray32", colSig = "red3", colLine = "#ff000080",
  log = "x", cex=0.45, xlab="mean expression", ylab="log fold change", ... )

\S4method{plotMA}{ExonCountSet}( object, FDR = 0.1, ... )

}
\arguments{
   \item{object}{
      either an ExonCountSet or a data.frame. If object is a data.frame, it must
      contain three columns, the first containing the mean expression values (for the
      x axis), the second the log fold change (for the y axis) and the third must be
      a logical vector indicating significance (for the coloring of the dots)
   }
   \item{FDR}{
      the false discovery rate, i.e., threshold to the adjusted p values, to be used 
      to colour the dots
   }
   \item{ylim}{
      The limits for the y axis. If missing, an attenpt is made to choose a sensible value.
      Dots exceeding the limits will be displayed as triangles at the limits, pointing outwards.
   }
   \item{colNonSig}{
      color to use for non-significant data points
   }
   \item{colSig}{
      color to use for significant data points      
   }
   \item{colLine}{
      color to use for the horizontal zero line
   }
   \item{log}{
      which axis should be logarithmic; will be passed to \code{\link{plot}}
   }
  \item{cex}{
     The \code{cex} parameter for \code{\link{plot}}.
   }
  \item{xlab}{
     The x axis label.
   }
  \item{ylab}{
     The y axis label.
   }
  \item{...}{
     Further parameters to be passed through to \code{\link{plot}}.
   }
}

\examples{
     \dontrun{
        data("pasillaExons", package="pasilla")
        pasillaExons <- estimateSizeFactors( pasillaExons )
        pasillaExons <- estimateDispersions( pasillaExons )
        pasillaExons <- fitDispersionFunction( pasillaExons )
        pasillaExons <- testForDEU( pasillaExons )
        pasillaExons <- estimatelog2FoldChanges( pasillaExons )
        plotMA( pasillaExons )
     }
}

