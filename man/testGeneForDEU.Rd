\name{testGeneForDEU}
\alias{testGeneForDEU}
\title{
      Test each a gene for Differential Exon Usage.
}
\description{
      For each exon of the gene, it applies a glm to a model frame of the gene using \code{modelFrameForGene}. It uses a negative binomial distribution for the glm using as a link function 1/dispersion, where dispersion is the highest value between the exon dispersion estimate and the common dispersion estimate. This function is designed to make a significance test between two fitted glm for each exon of a gene.
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
     Formula for the NULL model to fit a glm.  The values are going to be taken from the column names of the model frame generated from the function modelFrameForGene. As it is tested for each of the exons, a factor exonID can be added that will iterate over the exons of the gene fitting the glm. If it is left in NULL, the default formula is going to be "count~sample+exon+condition" as a NULL model.
  }
  \item{formula1}{
     Formula for the FULL model to fit a glm.  The values are going to be taken from the column names of the model frame generated from the function modelFrameForGene. As it is tested for each of the exons, a factor exonID can be added that will iterate over the exons of the gene fitting the glm. If it is left in NULL, the default formula is going to be "count~sample+exon+condition*I(exon==exonID)" as a NULL model. 
  }
  \item{glm.control}{
     A list of arguments to be passed to the glm.control.
   }
}
   \value{
	A vector of p values of each of the exons of the given geneID.
}
\examples{
	data("pasillaExons", package="pasilla")
	\dontrun{pasillaExons <- estimateSizeFactors( pasillaExons ) }
	\dontrun{pasillaExons <- estimateDispersions( pasillaExons )}
	\dontrun{testGeneForDEU(pasillaExons, "FBgn0085442")}
}
\seealso{
   \code{profileLogLikelihood}, \code{estimateExonDispersionsForModelFrame}, \code{estimateCommonDispersion}, \code{commonDispersion}, \code{testForDEU}
}
