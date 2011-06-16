\name{testGeneForDEU}
\alias{testGeneForDEU}
\title{
      Test each a gene for Differential Exon Usage.
}
\description{
For each exon of the gene, it applies a glm to a model frame (from the function \code{modelFrameForGene}) of the gene. It uses a negative binomial distribution for the glm using as a link function 1/dispersion, where dispersion is the highest value between the exon dispersion estimate and the fitted mean dispersion estimate. This function is designed to make a significance test between two fitted glm for each exon of a gene.
}
\usage{
testGeneForDEU(ecs, gene, formula0=NULL, formula1=NULL, glm.control=list(maxit = 100, epsilon = 3e-04))
}
\arguments{
  \item{ecs}{
     An ExonCountSet object.
  }
  \item{gene}{
     A gene identificator.
  }
  \item{formula0}{
     Formula for the NULL model to fit a glm.  The factors of the formula are going to be taken from the column names of the model frame generated from the function modelFrameForGene. If it is left in NULL, the default formula is going to be "count~sample+exon+condition" as a NULL model.
  }
  \item{formula1}{
     Formula for the FULL model to fit a glm.  The values are going to be taken from the column names of the model frame generated from the function modelFrameForGene. If it is left in NULL, the default formula is going to be "count~sample+exon+condition*I(exon==exonID)" as a FULL model.  Because a glm is fitted for each one of the exons, a factor "exonID" will iterate over each of the exons of the geneID, e.g. If a geneID contains exons E01, E02, E03,...,EN and formula=NULL, the function will fit N glms, the last part of the formula will change in the iterations as follows: I(exon==E01), I(exon==E02), I(exon==E03),...,I(exon==EN). 
  }
  \item{glm.control}{
     A list of arguments to be passed to the function glm.control.
   }
}
   \value{
	A data frame with columns "deviance", "df" (degrees of freedom) and pvalues from the test.
}
\examples{
	data("pasillaExons", package="pasilla")
	\dontrun{pasillaExons <- estimateSizeFactors( pasillaExons ) }
	\dontrun{pasillaExons <- estimateDispersions( pasillaExons )}
	\dontrun{testGeneForDEU(pasillaExons, "FBgn0085442")}
}
\seealso{
   \code{profileLogLikelihood}, \code{estimateExonDispersionsForModelFrame}, \code{testForDEU}
}
