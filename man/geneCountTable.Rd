\name{geneCountTable}
\alias{geneCountTable}
\title{
   Makes a count table for genes.
}
\description{
   This function returns a read count table of non-negative integers an ExonCountSet where each row is a gene and columns are samples.
}
\usage{
geneCountTable(ecs)
}
\arguments{
  \item{ecs}{
     An ExonCountSet object.
}
}
\seealso{
   \code{DESeq}
}   
\examples{
   	data("pasillaExons", package="pasilla")
	head(geneCountTable(pasillaExons))
}   

