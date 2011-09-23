\name{estimateDispersions}
\docType{methods}
\alias{estimateDispersions,ExonCountSet-method}
\alias{estimateDispersions}
\title{Estimate exon dispersions}
\description{
      This function estimates the Cox-Reid dispersion values using information
      across genes and exons. It stored this estimates in \code{fData(ecs)$dispersion_CR_est}.
}
\usage{
\S4method{estimateDispersions}{ExonCountSet}( cds, formula=count ~ sample + condition*exon, initialGuess=.01, nCores=1, minCount=10, maxExon=70)
}
\arguments{
  \item{cds}{
     An ExonCountSet object.
   }
  \item{formula}{
     Formula used in the glm to calculate the dispersion values.  The factors on the formula must be present in the design columns of the ExonCountSet object..
   }
  \item{initialGuess}{
	An arbitrary initial guess to start with the iterations.
   }
  \item{nCores}{
	Number of cores to be used to estimate the dispersions.  Multicore package must be loaded in order to split the job in several cores.
   }
  \item{minCount}{
        Minimum number of counts to consider an exon in the tests. This significantly increases the speed of the dispersion estimations and testing for differential exon usage. This is supported by the fact that small count exons have a low probability of being called significant, so it will not affect the results.
   }
  \item{maxExon}{
        Genes with more exons than this parameter will be discarded from the analysis. This is a speed issue. Currently, time of dispersion estimations and testing for differential exon usage significantly increases with number of exons.
   }
}
\value{
	An ExonCountSet object with dispersion \code{featureData(cds)$dispersion_CR_est}) parameters filled).
}
\examples{
	\dontrun{data("pasillaExons", package="pasilla")}
	\dontrun{pasillaExons <- estimateSizeFactors( pasillaExons )}
	\dontrun{pasillaExons <- estimateDispersions( pasillaExons )}
}
