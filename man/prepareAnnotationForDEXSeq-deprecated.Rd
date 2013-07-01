\name{prepareAnnotationForDEXSeq-deprecated}
\alias{prepareAnnotationForDEXSeq-deprecated}
\alias{prepareAnnotationForDEXSeq}
\title{
   Prepare annotation transcriptDb object for DEXSeq.
}
\description{
   WARNING: This function is deprecated and it has been replaced by the function disjointExons, from the package
   GenomicFeatures.
}
\usage{
prepareAnnotationForDEXSeq( transcriptDb, aggregateGenes=FALSE, includeTranscripts=TRUE )
}
\arguments{
  \item{transcriptDb}{
     An transcriptDb object.
  }
  \item{aggregateGenes}{
     Logical. Indicates whether two or more genes sharing an exon should be merged into an 'aggregate gene'. If 'no', the exons that can not be assiged to a single gene are ignored.
  }
  \item{includeTranscripts}{
     Logical. Indicates whether the transcript information of each exon should be added.
  }
}
\author{
     From code kindly provided by Mike Love.
}
\value{
     A GRanges object.
}
\examples{
\dontrun{
   library(GenomicFeatures)
   hse <- makeTranscriptDbFromBiomart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
   exonicParts <- prepareAnnotationForDEXSeq( hse )
}
}
