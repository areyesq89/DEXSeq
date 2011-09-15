\name{fitDispersions}
\alias{fitDispersions}
\title{Fit the mean-variance function.}
\description{
      This function fit all the Cox-Reid dispersion values on a mean-variance dependent function to estimate 
      alpha1 and alpha0 in alpha(mu) = (alpha1/mu + alpha0), where mu is equal to the mean.  After this, the maximum
      value of each exon between the Cox-Reid and the mean fitted value will be stored in fData$dispersion.
}
\usage{
fitDispersions(ecs)
}
\arguments{
  \item{ecs}{
     An ExonCountSet object.
  }
}
\value{
	An ExonCountSet object with information of the fit included, as well as the fData(ecs)$dispersion filled.
}
\examples{
	data("pasillaExons", package="pasilla")
	pasillaExons <- estimateSizeFactors( pasillaExons )
	\dontrun{pasillaExons <- estimateDispersion( pasillaExons, fitDispersions=FALSE )}
	\dontrun{pasillaExons <- fitDispersions( pasillaExons )}
}
