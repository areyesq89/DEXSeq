\name{geneIDs}
\alias{geneIDs}
\title{
   Accesor for the geneIDs in an ExonCountSet object.
}
\description{
   This function is an access for the gene identificators for each of the rows in the count table.
}
\usage{
counts(ecs)
}
\arguments{
  \item{ecs}{
     An ExonCountSet object.
}
}
\examples{
   ecs <- makeExampleExonCountSet()
   head( geneIDs( ecs ) )
}   
