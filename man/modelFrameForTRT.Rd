\name{modelFrameForTRT}
\alias{modelFrameForTRT}
\title{
   Makes the model frame of an exonCountSet object for the TRT method.
}
\description{
   Creates a data frame containing the model frame for a gene with the columns 
   sample, exon, size factors, their respective counts and the design annotation.
}
\usage{
modelFrameForTRT(ecs)
}
\arguments{
  \item{ecs}{
     An ExonCountSet object.
  }
}
\value{
     A data frame containing the model frame for a gene.
}
\examples{
data("pasillaExons", package="pasilla")
modelFrameForTRT(pasillaExons)
}

