\name{fitDispersions}
\alias{fitDispersions}
\title{Fit the mean-variance function.}
\description{
      This function fitts a mean-variance function over the Cox-Reid dispersion estimates.
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
	An ExonCountSet object with information of the fit added, and the dispersion estimate used for the test determined (maximum between the value of the fit and the individual exon dispersion estimate).
}
\examples{
	data("pasillaExons", package="pasilla")
	pasillaExons <- estimateSizeFactors( pasillaExons )
	\dontrun{pasillaExons <- estimateDispersion( pasillaExons, fitDispersions=FALSE )}
	\dontrun{pasillaExons <- fitDispersions( pasillaExons )}
}
