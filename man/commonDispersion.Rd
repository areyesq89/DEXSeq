\name{commonDispersion}
\alias{commonDispersion}
\title{
   Accesor for the common dispersion estimate of an ExonCountSet object.
}
\description{
   This function access the common dispersion field of an ExonCountSet object.
}
\usage{
commonDispersion(ecs)
}
\arguments{
  \item{ecs}{
     An ExonCountSet object.
}
}
\examples{
   ecs <- makeExampleExonCountSet()
   ecs <- estimateDispersions(ecs)
   commonDispersion(ecs)
############################
   ecs <- makeExampleExonCountSet()
   commonDispersion(ecs) <- 0.01
}   
