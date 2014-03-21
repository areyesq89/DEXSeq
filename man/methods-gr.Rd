\name{findOverlaps}
\docType{methods}
\alias{subsetByOverlaps}
\alias{subsetByOverlaps-methods}
\alias{subsetByOverlaps,DEXSeqResults,GenomicRanges-method}
\alias{findOverlaps}
\alias{findOverlaps-methods}
\alias{findOverlaps,DEXSeqResults,GenomicRanges-method}
\title{Methods from GRanges overlaps for the DEXSeq object}
\description{
  This function generates an MA plot.
}
\usage{

  \S4method{subsetByOverlaps}{DEXSeqResults,GenomicRanges}( query, subject, maxgap = 0L, 
    minoverlap = 1L, type = c("any", "start", "end", "within", "equal"), ignore.strand = FALSE )

  \S4method{findOverlaps}{DEXSeqResults,GenomicRanges}( query, subject, maxgap = 0L, 
    minoverlap = 1L, type = c("any", "start", "end", "within", "equal"), ignore.strand = FALSE )

}
\arguments{
  \item{query}{
     Either a \code{DEXSeqResults} object.
  }
  \item{subject}{
      A GRanges or GRangesList object.
  }
  \item{maxgap, minoverlap, type}{
    See \code{\link[IRanges]{findOverlaps}} in the IRanges package for
    a description of these arguments.
  }
  \item{ignore.strand}{
    See \code{\link[GenomicRanges]{findOverlaps}} in the GenomicRanges package for
    a description of this.
  }

}

\examples{
     \dontrun{
       
     }
}

