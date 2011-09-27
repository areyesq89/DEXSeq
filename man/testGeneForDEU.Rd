\name{testGeneForDEU}
\alias{testGeneForDEU}
\title{
      Test each a gene for differential exon usage.
}
\description{
      For each exon of the gene, it applies a glm using a model frame 
      (from the function \code{modelFrameForGene}) of the gene. It 
      uses a negative binomial distribution for the glm using as a 
      link function 1/dispersion, where dispersion is the highest 
      value between the exon dispersion estimate and the fitted mean 
      dispersion estimate. This function is designed to make a 
      significance test between two fitted glm for each exon of a 
      gene.
}
\usage{
testGeneForDEU(ecs, gene, formula0=NULL, formula1=NULL, glm.control=list(maxit = 100, epsilon = 3e-04))
}
\arguments{
  \item{ecs}{
     An ExonCountSet object.
  }
  \item{gene}{
     A gene ID to be tested for differential exon usage.
  }
  \item{formula0}{
     Formula for the NULL model to fit a the glm.  The factors must be 
     present in the design columns of the ExonCountSet object. As it is 
     tested for each of the exons, a factor exonID can be added to the 
     formula, so that it will iterate over the exons of the gene fitting 
     the glm for each of them. If it is left in NULL, the default formula is 
     "count~sample+exon+condition" for the NULL model.
  }
  \item{formula1}{
     Same as formula0, but for the test model.  If it is left in NULL, the 
     default formula will be "count~sample+exon+condition*I(exon==exonID)".  
     If added a factor "exonID", it will iterate over each of the exons of 
     the geneID, e.g. If a geneID contains exons E01, E02, E03,...,EN, 
     and the function is left in the default formula, the function will fit N glms, the 
     last part of the formula will change in the iterations as follows: 
     I(exon==E01), I(exon==E02), I(exon==E03),...,I(exon==EN).
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
