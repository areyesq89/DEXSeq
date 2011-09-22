\name{estimateDispersions}
\docType{methods}
\alias{estimateDispersions,ExonCountSet-method}
\alias{estimateDispersions}
\title{Estimate exon dispersions}
\description{
      This function estimates the Cox-Reid dispersion values for each exon of each gene 
      and stores the results it in fData.  Then it calls fitDispersions function to fit all 
      the dispersion values on a mean-variance dependent function to estimate 
      alpha1 and alpha0 in alpha(mu) = (alpha1/mu + alpha0), where mu is equal to the mean.  After this, the maximum
      between the Cox-Reid and the mean fitted value will be stored in fData$dispersion. This
      last value will be used for further testing.
}
\usage{
\S4method{estimateDispersions}{ExonCountSet}( cds, formula=count ~ sample + condition*exon, initialGuess=.01, nCores=1, minCount=10, maxExon=70, fitDispersions=TRUE)
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
  \item{fitDispersions}{
	If TRUE, fitDispersions will be called inside this function. This will estimate the parameters alpha0 and alpha1 in alpha(mu)=(alpha1/mu + alpha0), then calculate and store in fData(ecs)$dispersion the maximum value between the fitted value and the Cox-Reid estimate. 
   }
}
\value{
	An ExonCountSet object with dispersion (\code{featureData(cds)$dispersion} and \code{featureData(cds)$dispersion_CR_est}) parameters filled).
}
\examples{
	\dontrun{data("pasillaExons", package="pasilla")}
	\dontrun{pasillaExons <- estimateSizeFactors( pasillaExons )}
	\dontrun{pasillaExons <- estimateDispersions( pasillaExons )}
}
