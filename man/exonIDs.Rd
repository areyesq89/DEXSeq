\name{exonIDs}
\alias{exonIDs}
\alias{exonIDs<-}
\title{
   Accesor for the exonIDs in an ExonCountSet object.
}
\description{
   This function is an accesor for the exon identificators for each of the rows in the count table.
   Please make sure that the exon IDs are ordered alphanumerically.
}
\usage{
exonIDs(ecs)
exonIDs(ecs) <- value
}
\arguments{
  \item{ecs}{
     An ExonCountSet object.
   }
  \item{value}{
     An vector of gene identificators, one for each of the rows of the count data.
   }
}
\examples{
	library(DEXSeq)
   	data("pasillaExons", package="pasilla")
	geneIDs(pasillaExons)
}
