\name{plotDispEsts-methods}
\docType{methods}
\alias{plotDispEsts}
\alias{plotDispEsts-methods}
\alias{plotDispEsts,ExonCountSet-method}
\title{Diagnostic plot for dispersion estimates and dispersion fit}
\usage{
\S4method{plotDispEsts}{ExonCountSet}( object, ymin = NULL, cex = 0.45, 
}

\arguments{
  \item{object}{
     An ExonCountSet object.
   }
  \item{ymin}{
     The lower limit of the y axis. If NULL, a sensible value will be calculated.
   }
  \item{cex}{
     The \code{cex} parameter for \code{\link{plot}}.
   }
}
\value{
   None. A plot is displayed.
}
\examples{
   data("pasillaExons", package="pasilla")
   pasillaExons <- estimateSizeFactors( pasillaExons )
   pasillaExons <- estimateDispersions( pasillaExons )
   plotDispEsts( pasillaExons )
}