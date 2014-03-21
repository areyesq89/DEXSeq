\name{DEXSeq-deprecated}
\alias{DEXSeq-deprecated}
\alias{design,ExonCountSet}
\alias{design,ExonCountSet-deprecated}
\alias{newExonCountSet}
\alias{DEUresultTable}
\alias{subsetByGenes}
\alias{geneCountTable}
\alias{estimateExonDispersionsForModelFrame_BM}
\alias{estimateDispersions_BM}
\alias{testGeneForDEU_BM}
\alias{testForDEU_BM}
\alias{makeCompleteDEUAnalysis_BM}
\alias{read.HTSeqCounts}
\alias{countTableForGene}
\alias{fitDispersionFunction}
\alias{estimatelog2FoldChanges}
\alias{modelFrameForGene}
\alias{buildExonCountSet}
\alias{constructModelFrame}
\alias{getCountVector}
\alias{estimateExonDispersion}
\alias{testExonForDEU}
\alias{doCompleteDEUAnalysis}
\alias{newExonCountSet-deprecated}
\alias{DEUresultTable-deprecated}
\alias{subsetByGenes-deprecated}
\alias{geneCountTable-deprecated}
\alias{estimateExonDispersionsForModelFrame_BM-deprecated}
\alias{estimateDispersions_BM-deprecated}
\alias{testGeneForDEU_BM-deprecated}
\alias{testForDEU_BM-deprecated}
\alias{makeCompleteDEUAnalysis_BM-deprecated}
\alias{read.HTSeqCounts-deprecated}
\alias{countTableForGene-deprecated}
\alias{fitDispersionFunction-deprecated}
\alias{estimatelog2FoldChanges-deprecated}
\alias{modelFrameForGene-deprecated}
\alias{buildExonCountSet-deprecated}
\alias{constructModelFrame-deprecated}
\alias{getCountVector-deprecated}
\alias{estimateExonDispersion-deprecated}
\alias{testExonForDEU-deprecated}
\alias{doCompleteDEUAnalysis-deprecated}
\title{
  This functions are deprecated and will become defunct.
}
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
    \item{newExonCountSet: \code{DEXseqDataSet}}
    \item{DEUresultTable: \code{DEXseqResults}}
    \item{subsetByGenes: }
    \item{geneCountTable: }
    \item{estimateExonDispersionsForModelFrame_BM: estimateDispersions-DEXSeqDataSet}
    \item{estimateDispersions_BM: estimateDispersions-DEXSeqDataSet}
    \item{testGeneForDEU_BM: \code{testForDEU}}
    \item{testForDEU_BM: \code{testForDEU}}
    \item{doCompleteDEUAnalysis: \code{DEXSeq}}
    \item{makeCompleteDEUAnalysis_BM: \code{DEXSeq}}
    \item{read.HTSeqCounts: \code{DEXSeqDataSetFromHTSeq}}
    \item{countTableForGene: }
    \item{fitDispersionFunction: estimateDispersions-DEXSeqDataSet}
    \item{estimatelog2FoldChange: \code{estimateExonFoldChanges}}
    \item{modelFrameForGene: }
    \item{buildExonCountSet: \code{DEXSeqDataSetFromSE}}
    \item{constructModelFrame: }
    \item{getCountVector: }
    \item{estimateExonDispersion: estimateDispersions-DEXSeqDataSet}
    \item{testExonForDEU: \code{testForDEU}}
    \item{doCompleteDEUAnalysis: \code{DEXSeq}}
  }
}
