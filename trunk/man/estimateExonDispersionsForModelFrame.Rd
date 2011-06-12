\name{estimateExonDispersionsForModelFrame}
   \alias{estimateExonDispersionsForModelFrame}
   \title{Estimates exon dispersions}
   \description{
	This function calculates the individual exon dispersions given a model frame. It uses \code{profileLogLikelihood} to optimize the dispersion parameters for each of the exons on the glm.
	}
\usage{
      estimateExonDispersionsForModelFrame(modelFrame, formula=NULL, initialGuess=0.01)
}
\arguments{
     \item{modelFrame}{
        Model frame provided by the function \code{modelFrameForGene}.
     }
     \item{formula}{
     Formula for the initial glm used to calculate the dispersion values.  The values are going to be taken from the column names of the model frame generated from the function modelFrameForGene. If it is left in NULL, the default formula is going to be "count~sample+condition*exon".
     }
     \item{initialGuess}{
        An arbitraty initial guess to initiate the optimization.
     }
}
\value{
	A vector of exon dispersions.
}
\examples{
	ecs <- makeExampleExonCountSet()
	ecs <- estimateSizeFactors(ecs)
	estimateExonDispersionsForModelFrame(modelFrameForGene(ecs, "SG7861"))
}
\seealso{
   \code{profileLogLikelihood}, \code{estimateExonDispersionsForModelFrame}, \code{estimateDispersions}
}
