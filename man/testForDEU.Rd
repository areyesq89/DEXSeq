\name{testForDEU}
\alias{testForDEU}
\title{Test for Differential Exon Usage.}
\description{
      This function tests for differential exon usage for each of 
      the exons in the object. It stores the results in the fields 
      \code{fData(ecs)$pvalue} and \code{fData(ecs)$padjust}.
}
\usage{
testForDEU( ecs, 
   nCores = 1, 
   formula0 = ~ sample + exon, 
   formula1 = ~ sample + exon + condition : exon, 
   dispColumn="dispersion")
}
\arguments{
  \item{ecs}{
     An ExonCountSet object.
  }
  \item{formula0}{
     Formula for the null model to be used in the GLM fit. 
  }
  \item{formula1}{
     Formula for the full model to be used in the GLM fit. 
  }
  \item{nCores}{
     Number of CPUcores to be used to calculate the $p$-values.  
     The \code{parallel} package must be loaded to use more than 1 core.
  }
  \item{dispColumn}{
     Column name from the feature data frame from where to take the dispersions.
  }
}
\details{
   The terms in the formulas must be columns of \code{design(ecs)}.
}
\value{
     An ExonCountSet object with \code{fData(ecs)$pvalue} and 
     \code{fData(ecs)$padjust} data slots filled.
}
\examples{
     \dontrun{
        data("pasillaExons", package="pasilla")
        pasillaExons <- estimateSizeFactors( pasillaExons )
        pasillaExons <- estimateDispersions() pasillaExons )
        pasillaExons <- fitDispersionFunction( pasillaExons )
        pasillaExons <- testForDEU( pasillaExons )
     }
}
