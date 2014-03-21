\name{featureCounts}
\alias{featureCounts}
\alias{featureIDs}
\alias{featureIDs<-}
\alias{exonIDs}
\alias{exonIDs<-}
\alias{groupIDs}
\alias{groupIDs<-}
\alias{geneIDs}
\alias{geneIDs<-}
\alias{sampleAnnotation}
\title{Accessor functions for DEXSeqDataSet details}
\description{
  Accessor functions of the DEXSeqDataSet object.
}
\usage{

  featureCounts(object, normalized=FALSE)

  featureIDs( object )
  featureIDs( object ) <- value

  exonIDs( object )
  exonIDs( object ) <- value

  groupIDs( object )
  groupIDs( object ) <- value

  geneIDs( object )
  geneIDs( object ) <- value

  sampleAnnotation( object )
  
}
\arguments{
  \item{object}{
    An DEXSeqDataSet object.
  }
  \item{value}{
    A character vector to replace previous values.
  }
  \item{normalized}{
    Logical indicating whether or not to
    divide the counts by the size factors or normalization
    factors before returning (normalization factors always
    preempt size factors)
  }
}
\value{
  'featureCounts' access the counts per exonic region or
  feature region names. 'featureIDs' and 'exonIDs' are accessor
  functions for the exon bin or features identifiers.
  'groupIDs' and 'geneIDs' are accessor functions for
  the character vector grouping the features, for example exonIDs
  from the same gene are grouped together by having the
  same value in the geneIDs. 'sampleAnnotation' is an accessor
  for the information from the samples.
}
\examples{
  data(pasillaDEXSeqDataSet, package="pasilla")
  head( featureCounts( dxd ) )
  head( featureIDs(dxd))
  head( exonIDs(dxd))
}
