setClass( "ExonCountSet",
   contains = "eSet",
   representation = representation(
      designColumns = "character",
      dispFitCoefs = "numeric",
      formulas = "list",
      annotationFile = "character"
   ),
   prototype = prototype( new( "VersionedBiobase",
      versions = c( classVersion("eSet"), ExonCountSet = "1.0.3" ) ) )
)

setMethod( "plotDispEsts", signature(object="ExonCountSet"),
  function( object, ... ){
    .Deprecated("plotDispEsts", msg="The method plotDispEsts for ExonCountSet objects has been deprecated, use the method for the DEXSeqDataSet object")
  }
)

setMethod( "plotMA", signature( object="ExonCountSet" ),
  function( object, ... ){
    .Deprecated("plotMA", msg="The method plotMA for ExonCountSet objects has been deprecated, use the method for the DEXSeqDataSet or the DEXSeqResults object")
  }
)

setMethod("counts", signature(object="ExonCountSet"),
  function( object, ... ) {
    .Deprecated("counts")
  }
)

setMethod("sizeFactors",  signature(object="ExonCountSet"),
  function( object, ... ) {
    .Deprecated("counts")
  }
)

setMethod("design", signature(object="ExonCountSet"),
  function( object, ... ) {
    .Deprecated("")
  }
)

setReplaceMethod("design", signature(object="ExonCountSet"),
  function( object, ... ) {
    .Deprecated("")
  }
)

newExonCountSet <- function( ... ){
  .Deprecated("DEXSeqDataSet")
}


DEUresultTable <- function( ... )
{
  .Deprecated("DEXSeqResults")
}

subsetByGenes <- function( ... ) {
  .Deprecated("")
}

geneCountTable <- function( ... ) {
  .Deprecated("")
}


estimateExonDispersionsForModelFrame_BM <- function( ... )
{
  .Deprecated("estimateDispersions for DEXSeqDataSet")
}


estimateDispersions_BM <- function( ... )
{
  .Deprecated("estimateDispersions for DEXSeqDataSet")
}


testGeneForDEU_BM <- function ( ... )
{
  .Deprecated("testForDEU")
}


testForDEU_BM <- function( ... )
{
  .Deprecated("testForDEU")
}

makeCompleteDEUAnalysis_BM <- function( ... )
{
  .Deprecated("DEXSeq")
}

read.HTSeqCounts <- function( ... )
{
  .Deprecated("DEXSeqDataSetFromHTSeq")
}

setMethod("estimateSizeFactors", signature(object="ExonCountSet"),
   function( object, ... ){
     .Deprecated("estimateSizeFactors for DEXSeqDataSet")
   }
)

countTableForGene <- function( ... ) {
   .Deprecated("")
}

fitDispersionFunction <- function( ... )
{
  .Deprecated("estimateDispersions")
}

estimatelog2FoldChanges <- function( ... )
{
  .Deprecated("estimateExonFoldChanges")   
}

modelFrameForGene <- function( ... ) {
  .Deprecated("")
}

buildExonCountSet <- function( ... ){
  .Deprecated("DEXSeqDataSetFromSE")
}

constructModelFrame <- function( ... ){
  .Deprecated("")
}

getCountVector <- function( ... ) {
  .Deprecated("")
}

estimateExonDispersion <- function( ... ){
  .Deprecated("estimateDispersionsGeneEst", package="DESeq2")
}

testExonForDEU <- function( ... ){
  .Deprecated("testForDEU")
}

setMethod("estimateDispersions", signature(object="ExonCountSet"),
  function( object, ... ){
    .Defunct("estimateDispersions", msg="the methods for the signature ExonCountSet have been deprecated, please use DEXSeqDataSet object instead")
  } 
)

doCompleteDEUAnalysis <- function( ... )
{
  .Deprecated("DEXSeq")
}


########## DEFUNCT

prepareAnnotationForDEXSeq <- function( ... )
{
  .Defunct("disjointExons", "GenomicFeatures")
}

countReadsForDEXSeq <- function ( ... )
{
  .Defunct("summarizedOverlaps", "")
}
