\name{estimateDispersions}
\docType{methods}
\alias{estimateDispersions,ExonCountSet-method}
\alias{estimateDispersions}
\title{Estimate dispersions}
\description{
      This function estimates for each counting bin
      of the \code{ExonCountSet} \code{object} a dispersion value. It stores
      these values in \code{fData(ecs)$dispBeforeSharing}.
}
\usage{
\S4method{estimateDispersions}{ExonCountSet}( object,
         formula = ~ sample + exon + condition : exon,
         nCores = 1,
         minCount = 10)
}

\arguments{
  \item{object}{
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
likelihood method of McCarthy et al. (Nucl Acid Res., 2012, 40:4288), which they 
devised for the \code{edgeR} package.
}
\examples{
if(suppressWarnings(require("pasilla", quietly=TRUE, character.only=TRUE))){

   data("pasillaExons", package="pasilla")
   pasillaExons <- estimateSizeFactors( pasillaExons )
   pasillaExons <- estimateDispersions( pasillaExons )
   head( fData(pasillaExons)$dispBeforeSharing )

}}
