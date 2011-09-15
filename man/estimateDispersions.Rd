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
\S4method{estimateDispersions}{ExonCountSet}(cds, formula=count ~ sample + condition * exon, file=NULL, initialGuess=.01, n.cores=1, quiet=FALSE, fitDispersions=TRUE)
}
\arguments{
  \item{cds}{
     An ExonCountSet object.
   }
  \item{formula}{
     Formula used in the glm to calculate the dispersion values.  The factors on the formula must be present in the design columns of the ExonCountSet object..
   }
    \item{file}{
      A file to write the progress of the profile log likelihoods of the dispersions estimations.
   }
  \item{initialGuess}{
	An arbitrary initial guess to start with the iterations.
   }
  \item{n.cores}{
	Number of cores to be used to estimate the dispersions.  Multicore package must be loaded in order to split the job in several cores.
   }
  \item{quiet}{
	If TRUE, no progress report is shown.
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
