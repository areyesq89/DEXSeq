\docType{methods}
\name{ExonCountSet-deprecated,methods}
\alias{plotMA,ExonCountSet}
\alias{plotMA,ExonCountSet-deprecated}
\alias{plotDispEsts,ExonCountSet-method}
\alias{plotDispEsts,ExonCountSet-method,deprecated}
\alias{sizeFactors,ExonCountSet-method}
\alias{sizeFactors,ExonCountSet-method,deprecated}
\title{Deprecated methods for ExonCountSet object}
\description{
  The ExonCountSet object has been deprecated and substituted by
  the DEXSeqDataSet object.  Therefore, all the functions and methods
  for the ExonCountSet object will be deprecated and defunct at the 
  next releases. These functions have been subsituted by new functions 
  that are compatible with the newer versions of DEXSeq.
}
\details{
  The replacements of deprecated function are summarized in the
  following items.
  \itemize{
    \item{plotDispEsts,ExonCountSet: plotDispEsts,DESeqDataSet}
    \item{plotMA-ExonCountSet: plotMA-DESeqDataSet and plotMA-DESeqResults}
    \item{counts-ExonCountSet: counts-DESeqDataSet, featureCounts-DEXSeqDataSet, counts-DEXSeqResults}
    \item{sizeFactors-ExonCountSet: sizeFactors-DESeqDataSet}
    \item{design-ExonCountSet: sampleAnnotation-DEXSeqDataSet}
    \item{estimateSizeFactors-ExonCountSet: estimateSizeFactors-DEXSeqDataSet}
    \item{estimateDispersions-ExonCountSet: estimateDispersions-DEXSeqDataSet}
  }
}
