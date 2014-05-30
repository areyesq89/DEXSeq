
estimateSizeFactors.DEXSeqDataSet <- function(object, locfunc=median, geoMeans) {
  normFactors <- estimateSizeFactorsForMatrix(featureCounts(object), locfunc, geoMeans=geoMeans)
  sizeFactors(object) <- rep(normFactors, 2)
  maxExons <- length( unique( object@modelFrameBM$exon ) )
  object@modelFrameBM$sizeFactor <- rep(normFactors, each=maxExons)
  object
}
  
setMethod("estimateSizeFactors", signature(object="DEXSeqDataSet"),
          estimateSizeFactors.DEXSeqDataSet)


estimateDispersions.DEXSeqDataSet <- 
  function( object, fitType=c("parametric","local","mean"),
    maxit=100, quiet=FALSE, formula=design(object), BPPARAM=MulticoreParam(workers=1))
{
  if (is.null(sizeFactors(object)) & is.null(normalizationFactors(object))) {
    stop("first call estimateSizeFactors or provide a normalizationFactor matrix before estimateDispersions")
  }
  if (!is.null(dispersions(object))) {
    if (!quiet) message("you had estimated dispersions, replacing these")
    mcols(object) <- mcols(object)[,!(mcols(mcols(object))$type %in% c("intermediate","results"))]
  }
  stopifnot(length(maxit)==1)
  fitType <- match.arg(fitType, choices=c("parametric","local","mean"))

  splitParts <- sort(
    rep(seq_len(BPPARAM$workers), 
    length.out=nrow(object) ) )
  splitObject <- split( object, splitParts )

  modelMatrix <- rmDepCols(
    model.matrix(formula, as.data.frame(colData(object))))
  
  splitObject <- bplapply( splitObject, 
      function(x){
        estimateDispersionsGeneEst(x, 
          maxit=maxit, quiet=quiet, 
          modelMatrix = modelMatrix, 
          niter = 10)}, 
    BPPARAM=BPPARAM )

  mergeObject <- do.call( rbind, splitObject )
  mcols(object) <- mcols( mergeObject )
  assays(object) <- assays(mergeObject)

  mcols(object)$baseMean <- mcols(object)$exonBaseMean
#  library(genefilter)
  mcols(object)$baseVar <- mcols(object)$exonBaseVar
  mcols(object)$allZero <- unname( rowSums( featureCounts(object)) == 0 |
      rowSums(counts(object, normalized = TRUE)[, colData(object)$exon == "others"]) ==0 )

  object <- estimateDispersionsFit(object, fitType=fitType, quiet=quiet)
  splitObject <- split( object, splitParts )
  
  splitObject <- bplapply( splitObject, 
      function(x){
        estimateDispersionsMAP(x, 
          maxit=maxit, 
          quiet=quiet, 
          modelMatrix=modelMatrix)
      }, 
    BPPARAM=BPPARAM )

  mcols(object) <- mcols( do.call( rbind, splitObject ) )
  mcols(object)$baseMean <- unname( rowMeans( counts(object) ) )
  mcols(object)$baseVar <- unname( rowVars( counts(object) ) )
  mcols(object)$dispersion <- pmin( mcols(object)$dispersion, ncol(object) )

  object

}

setMethod( "estimateDispersions", signature(object="DEXSeqDataSet"),
          estimateDispersions.DEXSeqDataSet )

plotDispEsts.DEXSeqDataSet <- function(object, ymin,
       genecol = "black", fitcol = "red", finalcol = "dodgerblue",
       legend=TRUE, xlab, ylab, log = "xy", cex = 0.45, ...)
{
  stopifnot(is(object, "DEXSeqDataSet"))
  mcols(object)$baseMean <- 
    mcols(object)$exonBaseMean
  mcols(object)$baseVar <- 
    mcols(object)$exonBaseVar
  plotDispEsts( as(object, "DESeqDataSet"), genecol=genecol, fitcol=fitcol,
    finalcol=finalcol, legend=legend, log=log, cex=cex, ...)
}

setMethod( "plotDispEsts", signature(object="DEXSeqDataSet"),
  plotDispEsts.DEXSeqDataSet )


plotMA.DEXSeqDataSet <- function( object, alpha=0.1, ylim=c(-2, 2), foldChangeColumn=NULL, ...){
  stopifnot( is(object, "DEXSeqDataSet") )
  hasResults <- attr(object, "results")
  dexseqResults <- mcols(object)[,elementMetadata( mcols( object ) )$type == "DEXSeq results"]
  if( ncol( dexseqResults ) == 0 ){
    stop("first call estimateExonFoldChanges")
  }
  if( is.null(foldChangeColumn)){
    y <- dexseqResults[,grep("log2fold", colnames(dexseqResults))[1]]
  }else{
    y <- dexseqResults[,foldChangeColumn]
  }
  x <- rowMeans( featureCounts(object, normalized=TRUE) )
  df <- data.frame( x, y, results(object)$padj < 0.1 )
  ylim=c(-2, 2)
  plotMA(df, ylim, ... )
}

setMethod("plotMA", signature(object="DEXSeqDataSet"),
          plotMA.DEXSeqDataSet)

plotMA.DEXSeqResults <- function(object, alpha=0.1, ylim=c(-2, 2), foldChangeColumn=NULL, ...){
  stopifnot( is(object, "DEXSeqResults") )
  x <- rowMeans( counts(object, normalized=TRUE) )
  dexseqResults <- object[,which( elementMetadata( object )$type == "DEXSeq results" )]
  if( is.null(foldChangeColumn)){
    y <- dexseqResults[,grep("log2fold", colnames(dexseqResults))[1]]
  }else{
    y <- dexseqResults[,foldChangeColumn]
  }
  df <- data.frame( x, y, object$padj < 0.1 )
  ylim=c(-2, 2)
  plotMA(df, ylim, ... )
}

setMethod("plotMA", signature(object="DEXSeqResults"),
          plotMA.DEXSeqResults)
    
show.DEXSeqResults <- function(object){
  cat("\n")
  cat(elementMetadata(object)[colnames(object) == "pvalue","description"])
  cat("\n\n")
  show(DataFrame(object))
}

setMethod("show", signature(object="DEXSeqResults"),
          show.DEXSeqResults )

counts.DEXSeqResults <- function(object, normalized=FALSE){
    if( normalized ){
        t(t( object$countData ) / object@sampleData$sizeFactor)
    }else{
       object$countData
   }
}

setMethod( "counts", signature(object="DEXSeqResults"),
          counts.DEXSeqResults )


subsetByOverlaps.DEXSeqResults <- function( query, subject, maxgap = 0L, minoverlap = 1L,
         type = c("any", "start", "end", "within", "equal"),
         ignore.strand = FALSE ){
  stopifnot( is( query, "DEXSeqResults") )
  genomicData <- query$genomicData
  overlaps <- findOverlaps( query=genomicData, subject=subject, maxgap=maxgap,
    minoverlap=minoverlap, type=type, ignore.strand=ignore.strand )
  query[queryHits( overlaps ),]
}

setMethod("subsetByOverlaps", signature(query="DEXSeqResults", subject="GenomicRanges"),
          subsetByOverlaps.DEXSeqResults)

findOverlaps.DEXSeqResults <- function( query, subject, maxgap = 0L, minoverlap = 1L,
         type = c("any", "start", "end", "within", "equal"),
         ignore.strand = FALSE ){
  stopifnot( is( query, "DEXSeqResults") )
  genomicData <- query$genomicData
  overlaps <- findOverlaps( query=genomicData, subject=subject, maxgap=maxgap,
    minoverlap=minoverlap, type=type, ignore.strand=ignore.strand)
  overlaps
}

setMethod("findOverlaps", signature(query="DEXSeqResults", subject="GenomicRanges"),
          findOverlaps.DEXSeqResults)


