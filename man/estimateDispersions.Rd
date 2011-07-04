\name{estimateDispersions}
\docType{methods}
\alias{estimateDispersions,ExonCountSet-method}
\alias{estimateDispersions}
\title{Estimate exon dispersions}
\description{
      Makes for every gene a call to the internal functions estimateExonDispersionsForModelFrame, fits a mean dependent regression on them.
}
\usage{
\S4method{estimateDispersions}{ExonCountSet}(cds, formula=NULL, file=NULL, initialGuess=.01, quiet=FALSE)
}
\arguments{
  \item{cds}{
     An ExonCountSet object.
   }
  \item{formula}{
     Formula for the initial glm used to calculate the dispersion values.  The factors on the formula must be present as column names in the model frame generated from the function modelFrameForGene. If it is left in NULL, the default formula is going to be "count~sample+condition*exon".
   }
    \item{file}{
    A file where to write the iteration progress with dispersions and profile log likelihood values for each of the exons of the genes.
   }
  \item{initialGuess}{
	An arbitrary initial guess to start with the iterations.
   }
  \item{quiet}{
	If TRUE, no progress report is shown.
   }
}
\value{
	An ExonCountSet object with dispersion (\code{featureData(cds)$dispersion} and \code{featureData(cds)$dispersion_CR_est}) parameters filled as a featureData field for the Cox Reid dispersions as well as the dispersion values that will be used for the test (maximum between Cox-Reid dispersion and mean fitted dispersion).
}
\examples{
	\dontrun{data("pasillaExons", package="pasilla")}
	\dontrun{pasillaExons <- estimateSizeFactors( pasillaExons )}
	\dontrun{pasillaExons <- estimateDispersions( pasillaExons )}
}
