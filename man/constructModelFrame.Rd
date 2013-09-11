\name{constructModelFrame}
\alias{constructModelFrame}
\title{
   Returns the model frame of an exonCountSet object for the GLM fits.
}
\description{
   Creates a data frame containing the model frame for a gene with the columns 
   sample, exon, size factors, their respective counts and the design annotation.
}
\usage{
constructModelFrame(ecs)
}
\arguments{
  \item{ecs}{
     An ExonCountSet object.
  }
}
\value{
     A data frame containing the model frame for a gene.
}
\details{
	This function is called internally by several DEXSeq function, but is exposed as
	it might occasionally be usefule to users, too.
}
\examples{
data("pasillaExons", package="pasilla")
constructModelFrame(pasillaExons)
}

