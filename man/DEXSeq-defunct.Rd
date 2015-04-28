\name{DEXSeq-defunct}
\alias{DEXSeq-defunct}
\title{
  This functions are deprecated and will become defunct.
}
\description{
  The ExonCountSet object has been deprecated and substituted by
  the DEXSeqDataSet object.  Therefore, all the functions and methods
  for the ExonCountSet object were deprecated. These functions have 
  been subsituted by new functions and methods designed for DEXSeqDataSet objects.
}
\details{
  The replacements of defunct functions are summarized in the
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
