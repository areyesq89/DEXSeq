\name{estimateDispersionsTRT}
\alias{estimateDispersionsTRT}
\title{Estimate exon dispersions using the TRT method}
\description{
      This function estimates for each counting bin
      of the \code{ExonCountSet} \code{object} a dispersion value. It stores
      these values in \code{fData(ecs)$dispersionBeforeSharing}.
}
\usage{
estimateDispersionsTRT( ecs,
         formula= ~ sample + condition * exon,
         nCores=1,
         minCount=10)
}
\arguments{
  \item{ecs}{
     An ExonCountSet object.
   }
  \item{formula}{
     Formula used in the GLM to estimate the dispersion values.
     The terms in the formula must be design
     columns of the ExonCountSet object, the l.h.s. will be the counts
     for each exon.
   }
  \item{nCores}{
     Number of cores to be used to estimate the dispersions.
     The \code{parallel} package must be loaded in order to spread the
     job onto several cores.
   }
  \item{minCount}{
     Counting bins with less than \code{minCount} counts (summed over all
     samples) are skipped in the tests. This reduces computation
     time, as counting bins with very few counts cannot give
     a significant signal anyway. For skipped counting bins, the
     \code{testable} column in \code{fData} is set to \code{FALSE}.
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
   pasillaExons <- estimateDispersionsTRT( pasillaExons )

}}
