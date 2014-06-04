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
    .Defunct("plotDispEsts", msg="The method plotDispEsts for ExonCountSet objects has been deprecated, use the method for the DEXSeqDataSet object")
  }
)

setMethod( "plotMA", signature( object="ExonCountSet" ),
  function( object, ... ){
    .Defunct("plotMA", msg="The method plotMA for ExonCountSet objects has been deprecated, use the method for the DEXSeqDataSet or the DEXSeqResults object")
  }
)

setMethod("counts", signature(object="ExonCountSet"),
  function( object, ... ) {
    .Defunct("counts")
  }
)

setMethod("sizeFactors",  signature(object="ExonCountSet"),
  function( object, ... ) {
    .Defunct("counts")
  }
)

setMethod("design", signature(object="ExonCountSet"),
  function( object, ... ) {
    .Defunct("")
  }
)

setReplaceMethod("design", signature(object="ExonCountSet"),
  function( object, ... ) {
    .Defunct("")
  }
)

newExonCountSet <- function( ... ){
  .Defunct("DEXSeqDataSet")
}


DEUresultTable <- function( ... )
{
  .Defunct("DEXSeqResults")
}

subsetByGenes <- function( ... ) {
  .Defunct("")
}

geneCountTable <- function( ... ) {
  .Defunct("")
}


estimateExonDispersionsForModelFrame_BM <- function( ... )
{
  .Defunct("estimateDispersions for DEXSeqDataSet")
}


estimateDispersions_BM <- function( ... )
{
  .Defunct("estimateDispersions for DEXSeqDataSet")
}


testGeneForDEU_BM <- function ( ... )
{
  .Defunct("testForDEU")
}


testForDEU_BM <- function( ... )
{
  .Defunct("testForDEU")
}

makeCompleteDEUAnalysis_BM <- function( ... )
{
  .Defunct("DEXSeq")
}

read.HTSeqCounts <- function( ... )
{
  .Defunct("DEXSeqDataSetFromHTSeq")
}

setMethod("estimateSizeFactors", signature(object="ExonCountSet"),
   function( object, ... ){
     .Defunct("estimateSizeFactors for DEXSeqDataSet")
   }
)

countTableForGene <- function( ... ) {
   .Defunct("")
}

fitDispersionFunction <- function( ... )
{
  .Defunct("estimateDispersions")
}

estimatelog2FoldChanges <- function( ... )
{
  .Defunct("estimateExonFoldChanges")   
}

modelFrameForGene <- function( ... ) {
  .Defunct("")
}

buildExonCountSet <- function( ... ){
  .Defunct("DEXSeqDataSetFromSE")
}

constructModelFrame <- function( ... ){
  .Defunct("")
}

getCountVector <- function( ... ) {
  .Defunct("")
}

estimateExonDispersion <- function( ... ){
  .Defunct("estimateDispersionsGeneEst", package="DESeq2")
}

testExonForDEU <- function( ... ){
  .Defunct("testForDEU")
}

setMethod("estimateDispersions", signature(object="ExonCountSet"),
  function( object, ... ){
    .Defunct("estimateDispersions", msg="the methods for the signature ExonCountSet have been deprecated, please use DEXSeqDataSet object instead")
  } 
)

doCompleteDEUAnalysis <- function( ... )
{
  .Defunct("DEXSeq")
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
