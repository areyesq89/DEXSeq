\name{estimateExonDispersionsForModelFrame}
   \alias{estimateExonDispersionsForModelFrame}
   \title{Estimates exon dispersions}
   \description{
       This function calculates the individual profile log likelihood
       Cox-Reid exon dispersions for a single gene. It inputs the a 
       model frame generate from the function \code{\link{modelFrameForGene}}.
	}
\usage{
      estimateExonDispersionsForModelFrame(modelFrame, formula=NULL, mm=NULL, muhat=NULL, initialGuess=0.01)
}
\arguments{
     \item{modelFrame}{
        Model frame provided by the function \code{\link{modelFrameForGene}}.
     }
     \item{formula}{
        Formula for the glm used to calculate the Cox-Reid dispersions. The factors in the formula
        must be present in the column names of the model frame. If it is left in NULL, the default 
       formula will be "count~sample+condition*exon".
     }
     \item{mm}{
        A model matrix for the model frame, if NULL a model matrix will 
        be created from the parameters "formula" and "modelFrame".
     }
     \item{muhat}{
        Initial values for the coefficients in the optimization.  
        If NULL, initial values will be calculated using the "initialGuess".
     }
     \item{initialGuess}{
        An arbitrary initial guess to initiate the optimization.
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
        \code{\link{estimateDispersions}}
}
