\name{exonIDs}
\alias{exonIDs}
\title{
   Accesor for the exonIDs in an ExonCountSet object.
}
\description{
   This function is an accesor for the exon identificators for each of the rows in the count table. Please make sure that the exon identificators within a gene are ordered alphanumerically.
}
\usage{
exonID(ecs)
}
\arguments{
  \item{ecs}{
     An ExonCountSet object.
   }
}
\examples{
   ecs <- makeExampleExonCountSet()
   head( exonIDs( ecs ) )
}
