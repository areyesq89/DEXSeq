\name{modelFrameForGene}
\alias{modelFrameForGene}
\title{
   Makes the model frame for a geneID.
}
\description{
   Creates a data.frame containing the model frame for a gene with the columns, sample, exon, condition, sizeFactors and its respective count.
}
\usage{
modelFrameForGene(ecs, geneID)
}
\arguments{
  \item{ecs}{
     An ExonCountSet object.
  }
  \item{geneID}{
     A gene identificator contained in the ExonCountSet object.
  }
}
\value{
	A data.frame containing the model frame for a gene.
}
\examples{
   ecs <- makeExampleExonCountSet()
   modelFrameForGene(ecs, "SG10020")
}

