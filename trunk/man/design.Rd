\name{design}
\alias{design}
\title{
    Accessor function for the design information in a ExonCountSet object
}
\description{
    The design vector is a factor that assigns to each column of the count data
    a condition (or treatment, or phenotype, or the like). This information is stored
    in the ExonCountSet's "phenoData" slot as a row named "design".
}
\usage{
design(ecs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ecs}{
      An ExonCountSet
}
}
\author{
   Simon Anders, sanders@fs.tum.de
}
\examples{
   ecs <- makeExampleExonCountSet()
   design( ecs )
}
