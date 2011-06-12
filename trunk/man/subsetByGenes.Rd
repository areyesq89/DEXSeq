\name{subsetByGenes}
\alias{subsetByGenes}
\title{
    Making a subsetting the ExonCountSet object by geneIDs.
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
	ecs <- makeExampleExonCountSet()
	ecs2<-subsetBygenes(ecs, sample(levels(geneIDs(ecs)), 10))
}

