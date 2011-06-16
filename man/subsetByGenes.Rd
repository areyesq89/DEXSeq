\name{subsetByGenes}
\alias{subsetByGenes}
\title{
    Making an ExonCountSet out of a bigger ExonCountSet object given a subset of geneIDs.
}
\description{
    Generates an ExonCountSet object containing a subset of genes from another ExonCountSet.
}
\usage{
subsetByGenes(ecs, genes)
}
\arguments{
  \item{ecs}{
      An ExonCountSet
  }
  \item{genes}{
      Subset of geneIDs used to generate the subset ExonCountSet.
  }
}
\examples{
	data("pasillaExons", package="pasilla")
	ecs <- subsetByGenes(pasillaExons, sample(unique(geneIDs(pasillaExons)), 10))
}

