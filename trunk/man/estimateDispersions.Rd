\name{estimateDispersions}
\alias{estimateDispersions}
\title{Estimate exon dispersions and the common dispersion estimate.}
\description{
      Makes for every gene a call to the internal functions \code{estimateExonDispersionsForModelFrame} and \code{estimateCommonDispersion}.
}
\usage{
ecs <- estimateDispersions(ecs, formula=NULL, what="both")
}
\arguments{
  \item{ecs}{
     An ExonCountSet object.
  }
  \item{formula}{
     Formula for the initial glm used to calculate the dispersion values.  The values are going to be taken from the column names of the model frame generated from the function modelFrameForGene. If it is left in NULL, the default formula is going to be "count~sample+condition*exon".
  }
  \item{what}{
     The possible values are "both", "common", "exon".  Determines which of the dispersion parameters is to be calculated, it can me the common dispersion, the individual exon dispersions or both.
  }
}
\value{
	An ExonCountSet object with dispersion parameters filled as a featureData field for the exon dispersions and an environmental variable for the common dispersion estimate.
}
\seealso{
   \code{profileLogLikelihood}, \code{estimateExonDispersionsForModelFrame}, \code{estimateCommonDispersion}, \code{commonDispersion}, \code{featureData}
}
\examples{
	ecs <- makeExampleExonCountSet(ecs)
	ecs <- estimateSizeFactors(ecs)  
	ecs <- estimateDispersions(ecs)
	featureData(ecs)$dispersion
	commonDispersion(ecs)
}
