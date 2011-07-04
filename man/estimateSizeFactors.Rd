\name{estimateSizeFactors}
\docType{methods}
\alias{estimateSizeFactors,ExonCountSet-method}
\alias{estimateSizeFactors}
\title{
   Estimate the size factors for an ExonCountSet
}
\description{
   Estimate the size factors for an ExonCountSet
}
\usage{
\S4method{estimateSizeFactors}{ExonCountSet}(cds)
}
\arguments{
  \item{cds}{
     An ExonCountSet object
}
}
\details{
   This estimates the size factors and stores the information in the object. This function takes the count data from an ExonCountSet object (cds), estimates the size factors 
   as follows: Each column (sample) is divided by the geometric means of the rows. The median
   of these ratios (skipping the genes with a geometric mean of zero) is used as the
   size factor for this column.
}
\value{
   The ExonCountSet passed as parameters, with the size factors filled in.
}
\author{
   Simon Anders, sanders@fs.tum.de
}
\examples{
	data("pasillaExons", package="pasilla")
	pasillaExons <- estimateSizeFactors( pasillaExons )
        sizeFactors( pasillaExons )
}
