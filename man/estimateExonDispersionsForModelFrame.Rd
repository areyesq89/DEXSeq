\name{estimateExonDispersionsForModelFrame}
   \alias{estimateExonDispersionsForModelFrame}
   \title{Estimates exon dispersions}
   \description{
	This function calculates the individual exon dispersions given a model frame. It uses \code{profileLogLikelihood} to optimize the dispersion parameters for each of the exons on the glm.
	}
\usage{
      estimateExonDispersionsForModelFrame(modelFrame, formula=NULL, mm=NULL, muhat=NULL, initialGuess=0.01)
}
\arguments{
     \item{modelFrame}{
        Model frame provided by the function \code{modelFrameForGene}.
     }
     \item{formula}{
     Formula for the initial glm used to calculate the dispersion values.  The values are going to be taken from the column names of the model frame generated from the function modelFrameForGene. If it is left in NULL, the default formula is going to be "count~sample+condition*exon".
     }
     \item{mm}{
      A model matrix for the model frame, if NULL a model matrix will be created with the options "formula" and "modelFrame".
     }
     \item{muhat}{
      Initial values for the coefficients in the optimization.  If NULL, initial values will be calculated using the "initialGuess" dispersion value.
     }
     \item{initialGuess}{
        An arbitraty initial guess to initiate the optimization.
     }
}
\value{
	A vector of exon dispersions.
}
\examples{
	data("pasillaExons", package="pasilla")
	pasillaExons <- estimateSizeFactors( pasillaExons )
	estimateExonDispersionsForModelFrame(modelFrameForGene(pasillaExons, "FBgn0085442"))
}
\seealso{
   \code{profileLogLikelihood}, \code{estimateExonDispersionsForModelFrame}, \code{estimateDispersions}
}
