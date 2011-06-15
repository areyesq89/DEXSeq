\name{estimateDispersions}
\alias{estimateDispersions}
\title{Estimate exon dispersions}
\description{
      Makes for every gene a call to the internal functions estimateExonDispersionsForModelFrame, fits a mean dependent regression on them.
}
\usage{
estimateDispersions(ecs, formula=NULL, file=NULL, initialGuess=.01, quiet=FALSE)
}
\arguments{
  \item{ecs}{
     An ExonCountSet object.
   }
  \item{formula}{
     Formula for the initial glm used to calculate the dispersion values.  The values are going to be taken from the column names of the model frame generated from the function modelFrameForGene. If it is left in NULL, the default formula is going to be "count~sample+condition*exon".
   }
    \item{file}{
    A file where to write the iterations and dispersion and profile log likelihood values foreach of the exons of the genes.
   }
  \item{initialGuess}{
	An arbitrary initial guess to start with the iterations.
   }
  \item{quiet}{
	If TRUE, no progress report is shown.
   }
}
\value{
	An ExonCountSet object with dispersion parameters filled as a featureData field for the Cox Reid dispersions as well as the fitted dispersions.
}
\examples{
	\dontrun{data("pasillaExons", package="pasilla")}
	\dontrun{pasillaExons <- estimateSizeFactors( pasillaExons )}
	\dontrun{pasillaExons <- estimateDispersions( pasillaExons )}
}
