\name{testForDEU}
\alias{testForDEU}
\title{Test for Differential Exon Usage.}
\description{
      This function will test for differential exon usage for each of the genes in the object, it will store the results in the fields of \code{fData(ecs)$pvalue} and \code{fData(ecs)$padjust}.
}
\usage{
testForDEU(ecs, formula0=NULL, formula1=NULL, padjust=TRUE, n.cores=1, quiet=FALSE)
}
\arguments{
  \item{ecs}{
     An ExonCountSet object.
  }
  \item{formula0}{
     Formula for the NULL model to fit a glm.  The factors must be present in the design columns of the ExonCountSet object. As it is tested for each of the exons, a factor exonID can be added to the formula, so that it will iterate over the exons of the gene fitting the glm. If it is left in NULL, the default formula is "count~sample+exon+condition" for the NULL model.
  }
  \item{formula1}{
     Same as formula0, but for the test model.  If it is left in NULL, the default formula will be "count~sample+exon+condition*I(exon==exonID)".  A factor "exonID" will iterate over each of the exons of the geneID, e.g. If a geneID contains exons E01, E02, E03,...,EN and formula=NULL, the function will fit N glms, the last part of the formula will change in the iterations as follows: I(exon==E01), I(exon==E02), I(exon==E03),...,I(exon==EN).
  }
  \item{padjust}{
     If TRUE, multiple testing correction will be included in the field featureData(ecs)$padjust.
  }
  \item{n.cores}{
	Number of cores to be used to estimate the dispersions.  multicore package must be loaded in order to split the job in several cores.
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
