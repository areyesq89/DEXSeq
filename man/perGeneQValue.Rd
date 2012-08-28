\name{perGeneQValue}
\alias{perGeneQValueBySimulation}
\alias{perGeneQValueExact}
\alias{perGeneQValue}
\title{Summarize per-exon p-values into per-gene q-values.}
\description{
  The use case for this function is the following analysis: given
  per-exon p-values for null hypothesis H0, we can determine the number
  of genes in which at least for one exon H0 is rejected. What is the
  associated false disovery rate?.
}
\usage{
perGeneQValue(ecs, p = "pvalue", method = perGeneQValueExact)
perGeneQValueBySimulation(pGene, theta, geneSplit, nperm = 24)
perGeneQValueExact(pGene, theta, geneSplit)
}
\arguments{
\item{ecs}{An \code{ExonCountSet} object. \code{fData(ecs)} is required
  to have columns \code{testable} and \code{geneID}.}
\item{p}{A character string indicating the name of the slot in
  \code{fData(ecs)} from which to take the per-exon p-values.}
}
\details{Details
}
\value{
  A named numeric vector, values are per-gene q-values, names are gene.
}
\examples{
 ## example code
}
\seealso{See also
}
