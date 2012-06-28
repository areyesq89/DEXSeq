\name{estimateDispersions}
\docType{methods}
\alias{estimateDispersions,ExonCountSet-method}
\alias{estimateDispersions}
\title{Estimate exon dispersions}
\description{
      This function estimates for each counting bin
      of the \code{ExonCountSet} \code{object} a dispersion value. It stores
      these values in \code{fData(ecs)$dispersionBeforeSharing}.
}
\usage{
\S4method{estimateDispersions}{ExonCountSet}( object,
         formula=count ~ sample + condition * exon,
         initialGuess=.01, nCores=1, minCount=10,
         maxExon=70, quiet=FALSE, file="")
}
\arguments{
  \item{object}{
     An ExonCountSet object.
   }
  \item{formula}{
     Formula used in the GLM to estimate the dispersion values.
     The terms in the formula must be design
     columns of the ExonCountSet object, the l.h.s. must be \code{count}.
   }
  \item{initialGuess}{
     An initial guess for the dispersion values to initiate the optimization.
   }
  \item{nCores}{
     Number of cores to be used to estimate the dispersions.
     The \code{multicore} package must be loaded in order to spread the
     job onto several cores.
   }
  \item{minCount}{
     Counting bins with less than \code{minCount} counts (summed over all
     samples) are skipped in the tests. This reduces computation
     time, as counting bins with very few counts cannot give
     a significant signal anyway. For skipped counting bins, the
     \code{testable} column in \code{fData} is set to \code{FALSE}.
   }
  \item{maxExon}{
     Genes with more than \code{maxExon} counting bins will be skipped in the
     test. This option can be useful when otherwise genes with very many counting bins
     use up extremely long computation time for dispersion estimation and testing for
     differential exon usage.
   }
   \item{quiet}{
     If \code{TRUE}, no progress report is shown.  In case the session is
     not an interactive session and a progress report is wanted, include
     a file name in the parameter \code{file}.
   }
   \item{file}{
     A file name to write the progress reports. If \code{file} is
     \code{""}, the output
     will be written to the standard output connection.
   }
}
\value{
     An object of class \code{ExonCountSet} with dispersion \code{featureData(object)$dispersion_CR_est})
     parameters filled).
}
\details{For the dispersion estimation, we use the Cox-Reid conditional maximum
likelihood method of Gordon Smyth et al., which they devised for the \code{edgeR}
package.
}
\examples{
if(suppressWarnings(require("pasilla", quietly=TRUE, character.only=TRUE))){

   data("pasillaExons", package="pasilla")
   pasillaExons <- estimateSizeFactors( pasillaExons )
   pasillaExons <- estimateDispersions( pasillaExons )

}}
