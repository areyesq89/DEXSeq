\name{testForDEU}
\alias{testForDEU}
\title{Test for Differential Exon Usage.}
\description{
      This function makes a call to \code{testGeneForDEU} for each of the geneIDs, stores in the fields of \code{fData(ecs)$pvalue}, \code{fData(ecs)$pll} and \code{fData(ecs)$df}
}
\usage{
testForDEU(ecs, formula0=NULL, formula1=NULL, padjust=TRUE, quiet=FALSE)
}
\arguments{
  \item{ecs}{
     An ExonCountSet object.
  }
  \item{formula0}{
     Formula for the NULL model to fit a glm.  The values are going to be taken from the column names of the model frame generated from the function modelFrameForGene. As it is tested for each of the exons, a factor exonID can be added that will iterate over the exons of the gene fitting the glm. If it is left in NULL, the default formula is "count~sample+exon+condition" for the NULL model.
  }
  \item{formula1}{
     Formula for the FULL model to fit a glm.  The values are going to be taken from the column names of the model frame generated from the function modelFrameForGene. If it is left in NULL, the default formula is going to be "count~sample+exon+condition*I(exon==exonID)" as a FULL model.  Because a glm is fitted for each one of the exons, a factor "exonID" will iterate over each of the exons of the geneID, e.g. If a geneID contains exons E01, E02, E03,...,EN and formula=NULL, the function will fit N glms, the last part of the formula will change in the iterations as follows: I(exon==E01), I(exon==E02), I(exon==E03),...,I(exon==EN). 
  }
  \item{padjust}{
     Logical argument indicating whether to include the multiple testing correction in the field featureData(ecs)$padjust.
  }
  \item{quiet}{
    If TRUE, no progress report is shown. 
  }
}
\value{
	An ExonCountSet object with \code{fData(ecs)$pvalue} and \code{fData(ecs)$padjust} data slots filled.
}
\examples{
	data("pasillaExons", package="pasilla")
	pasillaExons <- estimateSizeFactors( pasillaExons )
	\dontrun{pasillaExons <- estimateDispersion( pasillaExons )}
	\dontrun{pasillaExons <- testForDEU( pasillaExons )}
}
\seealso{
   \code{profileLogLikelihood}, \code{estimateExonDispersionsForModelFrame}
}
