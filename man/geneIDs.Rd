\name{geneIDs}
\alias{geneIDs}
\alias{geneIDs<-}
\title{
   Accesor for the geneIDs in an ExonCountSet object.
}
\description{
   This function is an access for the gene identificators for each of the rows in the count table.
}
\usage{
geneIDs(ecs)
geneIDs(ecs) <- value
}
\arguments{
  \item{ecs}{
     An ExonCountSet object.
  }
  \item{value}{
     An vector of exon identificators, one for each of the rows of the count data.
   }
}
\examples{
library(DEXSeq)
data("pasillaExons", package="pasilla")
head( geneIDs( pasillaExons ) )
}   
