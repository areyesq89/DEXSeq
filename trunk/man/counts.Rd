\name{counts}
\alias{counts}
\title{
   Accessor for the 'counts' slot of a ExonCountSet object.
}
\description{
   The counts slot holds the count data as a matrix of non-negative integer
   count values, one row for each observational unit (an exon of a gene), and
   one column for each sample.
}
\usage{
counts(ecs)
}
\arguments{
  \item{ecs}{
     An ExonCountSet object
}
}
\examples{
   ecs <- makeExampleCountDataSet()
   head( counts( ecs ) )
}   
