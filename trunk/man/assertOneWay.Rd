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
   ecs <- makeExampleExonCountSet()
   assertOneWay(ecs)
}

