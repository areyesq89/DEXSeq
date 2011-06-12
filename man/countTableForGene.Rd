\name{countTableForGene}
\alias{countTableForGene}
\title{
    Count table for a given geneID.
}
\description{
	This function returns a matrix of non negative integers containing a count table for a specified geneID from an ecs object.  The count table contains one row for every exon of the gene and a column for every sample.
}
\usage{
countTableForGene(ecs, geneID, normalized=FALSE)
}
\arguments{
  \item{ecs}{
      An ExonCountSet.
}
  \item{geneID}{
      A geneID to get the count table.
}
  \item{normalized}{
      Logic argument indicating whether to get the raw counts or the normalized counts by the size factor.
  }
}

\seealso{
\code{\link{estimateSizeFactors}}
}
\examples{
   ecs <- makeExampleExonCountSet()
   ecs <- estimateSizeFactors( ecs )
   countTableForGene(ecs, "SG7861", normalized=FALSE)
}


