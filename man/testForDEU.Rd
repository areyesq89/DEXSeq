\name{testForDEU}
\alias{testForDEU}
\title{Test for Differential Exon Usage.}
\description{
      This function makes a call to \code{testGeneForDEU} foreach of the geneIDs, stores in the fields of \code{fData(ecs)$pvalue}, \code{fData(ecs)$pll} and \code{fData(ecs)$df}
}
\usage{
testForDEU(ecs)
}
\arguments{
  \item{ecs}{
     An ExonCountSet object.
  }
  \item{padjust}{
     Logical argument indicating whether to include the multiple testing correction in the field featureData(ecs)$padjust.
  }
  \item{formula0}{
     Formula for the NULL model to fit a glm.  The values are going to be taken from the column names of the model frame generated from the function modelFrameForGene. As it is tested for each of the exons, a factor exonID can be added that will iterate over the exons of the gene fitting the glm. If it is left in NULL, the default formula is "count~sample+exon+condition" for the NULL model.
  }
  \item{formula1}{
     Formula for the FULL model to fit a glm.  The values are going to be taken from the column names of the model frame generated from the function modelFrameForGene. As it is tested for each of the exons, a factor exonID can be added that will iterate over the exons of the gene fitting the glm. If it is left in NULL, the default formula is "count~sample+exon+condition*I(exon==exonID)" for the FULL model. 
  }
}
\value{
	An ExonCountSet object with \code{fData(ecs)$pvalue} field filled with pvalues from testing each of the genes.
}
\examples{
	ecs <- makeExampleExonCountSet(ecs)
	ecs <- estimateSizeFactors(ecs)
	ecs <- estimateDispersion(ecs)
	ecs <- testForDEU(ecs)
}
\seealso{
   profileLogLikelihood, \code{estimateExonDispersionsForModelFrame}, \code{estimateCommonDispersion}, \code{commonDispersion}
}
