\name{assertOneWay}
\alias{assertOneWay}
\title{
   Assert one way.
}
\description{
   This function makes sure of the correct structure of the object for proper functioning of further functions.
}
\usage{
assertOneWay(ecs)
}
\arguments{
  \item{ecs}{
     An ExonCountSet object.
}
}
\examples{
   	data("pasillaExons", package="pasilla")
	design( pasillaExons) <- design(pasillaExons)[,1]
	assertOneWay(pasillaExons)
}

