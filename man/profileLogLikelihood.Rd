\name{profileLogLikelihood}
\alias{profileLogLikelihood}
\title{
   Estimates profile loglikelihood for the dispersion parameters.
}
\description{
   It is an internal use function that calculates the profile loglikelihood for the dispersion estimation parameters of the count values using the Cox-Reid estimator.
}
\usage{
profileLogLikelihood(disp, mm, y, muhat)
}
\arguments{
  \item{disp}{
     A dispersion value or vector of dispersion values to calculate the profile loglikelihood from.
  }
  \item{mm}{
     A model matrix provided by the function \code{model.matrix}.
  }
  \item{y}{
     A vector of count values to which the profile loglikelihood of its dispersion is being calculated.
  }
  \item{muhat}{
    A vector of values with the residuals of a previously called glm.
  }
}  
\seealso{
   \code{estimateCommonDispesion}, \code{estimateExonDispersionsForModelFrame}, \code{estimateDispersions}
}
