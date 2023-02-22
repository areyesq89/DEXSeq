
estimateSizeFactors.DEXSeqDataSet <- function(object, locfunc=median, geoMeans) {
  # Temporary hack for backward compatibility with "old" DEXSeqDataSet
  # objects. Remove once all serialized DEXSeqDataSet objects around have
  # been updated.
  if (!.hasSlot(object, "rowRanges"))
      object <- updateObject(object)
  normFactors <- estimateSizeFactorsForMatrix(featureCounts(object), locfunc, geoMeans=geoMeans)
  sizeFactors(object) <- rep(normFactors, 2)
  maxExons <- length( unique( object@modelFrameBM$exon ) )
  object@modelFrameBM$sizeFactor <- rep(normFactors, each=maxExons)
  object
}

setMethod("estimateSizeFactors", signature(object="DEXSeqDataSet"),
          estimateSizeFactors.DEXSeqDataSet)

estimateDispersions.DEXSeqDataSet <- function( object, fitType=c("parametric","local","mean", "glmGamPoi"),
      maxit=100, niter=10, quiet=FALSE, formula=design(object), BPPARAM=SerialParam()) {
      ## Temporary hack for backward compatibility with "old" DEXSeqDataSet
      ## objects. Remove once all serialized DEXSeqDataSet objects around have
      ## been updated.
    if (!.hasSlot(object, "rowRanges"))
        object <- updateObject(object)
    if (is.null(sizeFactors(object)) & is.null(normalizationFactors(object))) {
        stop("first call estimateSizeFactors or provide a normalizationFactor matrix before estimateDispersions")
    }
    if (!is.null(dispersions(object))){
        if (!quiet) message("you had estimated dispersions, replacing these")
        mcols(object) <- mcols(object)[,!(mcols(mcols(object))$type %in% c("intermediate","results"))]
    }
    stopifnot(length(maxit)==1)
    fitType <- match.arg(fitType, choices=c("parametric","local","mean", "glmGamPoi"))
    
    allVars <- all.vars(formula)
    if( any(!allVars %in% colnames( colData(object) )) ){
        notPresent <- allVars[!allVars %in% colnames( colData(object) )]
        notPresent <- paste(notPresent, collapse=",")
        stop(sprintf("the variables '%s' of the parameter 'formula' are not specified in the columns of the colData", notPresent ))
    }

    if( fitType == "glmGamPoi" & !is(BPPARAM, "SerialParam") ){
        message("\nParallelization has not been implemented for estimation of dispersions using glmGamPoi, using a single core for this estimation.\n")
        BPPARAM <- SerialParam()
    }

    if( is( BPPARAM, "SerialParam" ) ){
        numParts <- 1L
    }else{
        numParts <- BPPARAM$workers
    }
    
    splitParts <- sort( rep( seq_len( numParts ), length.out=nrow(object) ) )
    splitObject <- split( object, splitParts )
    
    modelMatrix <- rmDepCols(
        model.matrix(formula, as.data.frame(colData(object))))
    
    glmType <- ifelse( fitType == "glmGamPoi", "glmGamPoi", "DESeq2" )
    
    splitObject <- bplapply( splitObject,
        function(x, ... ){
            estimateDispersionsGeneEst(x,
                maxit=maxit, quiet=quiet,
                modelMatrix = modelMatrix,
                niter = niter, type=glmType ) },
        maxit=maxit, quiet=quiet,
        modelMatrix=modelMatrix,
        glmType=glmType,
        niter=niter,
        BPPARAM=BPPARAM )

    mergeObject <- do.call( rbind, splitObject )
    matchedNames <- match( rownames(object), rownames(mergeObject))
    mcols(object) <- mcols( mergeObject )[matchedNames,]
    assays(object) <- assays(mergeObject[matchedNames,])
    
    mcols(object)$baseMean <- rowMeans( featureCounts(object, normalized=TRUE) )
    mcols(object)$baseVar <- mcols(object)$exonBaseVar
    mcols(object)$allZero <-
                    unname( rowSums( featureCounts(object)) == 0 |
                        rowSums(counts(object)[, colData(object)$exon == "others"]) ==0 )

    ## object <- estimateDispersionsFit(object, fitType=fitType, quiet=quiet)
    object <- estimateDispersionsFit(object, fitType=fitType, quiet=quiet)

    dispPriorVar <- estimateDispersionsPriorVar(object, modelMatrix=modelMatrix)

    splitObject <- split( object, splitParts )
    
    splitObject <- bplapply( splitObject,
        function(x, ... ){
            estimateDispersionsMAP(x,
                maxit=maxit,
                quiet=quiet,
                modelMatrix=modelMatrix,
                dispPriorVar=dispPriorVar,
                type=glmType)
        },
        maxit=maxit, quiet=quiet,
        modelMatrix=modelMatrix,
        dispPriorVar=dispPriorVar,
        glmType=glmType,
        BPPARAM=BPPARAM )

    if( fitType == "glmGamPoi" ){
        attr(object, "quasiLikelihood_df0") <- attr(splitObject[[1L]], "quasiLikelihood_df0")
    }
    
    mergeObject <- do.call( rbind, splitObject )
    matchedNames <- match( rownames(object), rownames(mergeObject) )
    mcols(object) <- mcols( mergeObject )[matchedNames,]
    mcols(object)$baseMean <- unname( rowMeans( counts(object, normalized=TRUE) ) )
    mcols(object)$baseVar <- unname( rowVars( counts(object, normalized=TRUE) ) )
    mcols(object)$dispersion <- pmin( mcols(object)$dispersion, ncol(object) )
    object
}

#' @docType methods
#' @name estimateDispersions
#' @aliases estimateDispersions,DEXSeqDataSet
#' @title Estimate the dispersions for a DEXSeqDataSet
#'
#' @param object A DEXSeqDataSet.
#' @param fitType Either "parametric", "local", "mean" or "glmGamPoi"
#' for the type of fitting of dispersions to the mean
#' intensity. See ?estimateDispersions,DESeqDataSet-method for details.
#' If "glmGamPoi" is selected, the GLM fitter from "glmGamPoi" is also
#' used for estimating the dispersion estimates.
#' @param  maxit Control parameter: maximum number of iterations to allow for convergence
#' @param niter Number of times to iterate between estimation of means and estimation of dispersion.
#' @param quiet Whether to print messages at each step.
#' @param formula Formula used to fit the dispersion estimates.
#' @param BPPARAM A "BiocParallelParam" instance. See \code{?bplapply} for details.
#'
#' @return A DEXSeqDataSet with the dispersion information filled in as metadata columns.
#' @description This function obtains dispersion estimates for negative binomial distributed data for the specific case for DEXSeq.
#'
#' @details See ?estimateDispersions,DESeqDataSet-method for details.
#'
#' @examples
#' data(pasillaDEXSeqDataSet, package="pasilla")
#' dxd <- estimateSizeFactors( dxd )
#' dxd <- estimateDispersions( dxd )
#'
#' @exportMethod estimateDispersions
#' 
setMethod( "estimateDispersions", signature(object="DEXSeqDataSet"),
          estimateDispersions.DEXSeqDataSet )

plotDispEsts.DEXSeqDataSet <- function(object, ymin,
       genecol = "black", fitcol = "red", finalcol = "dodgerblue",
       legend=TRUE, xlab, ylab, log = "xy", cex = 0.45, ...)
{
  stopifnot(is(object, "DEXSeqDataSet"))
  # Temporary hack for backward compatibility with "old" DEXSeqDataSet
  # objects. Remove once all serialized DEXSeqDataSet objects around have
  # been updated.
  if (!.hasSlot(object, "rowRanges"))
      object <- updateObject(object)
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
  # Temporary hack for backward compatibility with "old" DEXSeqDataSet
  # objects. Remove once all serialized DEXSeqDataSet objects around have
  # been updated.
  if (!.hasSlot(object, "rowRanges"))
      object <- updateObject(object)
  hasResults <- attr(object, "results")
  dexseqResults <- mcols(object)[,mcols( mcols( object ) )$type == "DEXSeq results"]
  if( ncol( dexseqResults ) == 0 ){
    stop("first call estimateExonFoldChanges")
  }
  if( is.null(foldChangeColumn)){
    y <- dexseqResults[,grep("log2fold", colnames(dexseqResults))[1]]
  }else{
    y <- dexseqResults[,foldChangeColumn]
  }
  x <- rowMeans( featureCounts(object, normalized=TRUE) )
  df <- data.frame( x, y, results(object)$padj < alpha )
#  ylim=c(-2, 2)
  plotMA(df, ylim, ... )
}

setMethod("plotMA", signature(object="DEXSeqDataSet"),
          plotMA.DEXSeqDataSet)

plotMA.DEXSeqResults <- function(object, alpha=0.1, ylim=c(-2,2), foldChangeColumn=NULL, ...){
  stopifnot( is(object, "DEXSeqResults") )
  # Temporary hack for backward compatibility with "old" DEXSeqDataSet
  # objects. Remove once all serialized DEXSeqDataSet objects around have
  # been updated.
  if (!.hasSlot(object, "rowRanges"))
      object <- updateObject(object)
  x <- rowMeans( counts(object, normalized=TRUE) )
  dexseqResults <- object[,which( mcols( object )$type == "DEXSeq results" )]
  if( is.null(foldChangeColumn)){
    y <- dexseqResults[,grep("log2fold", colnames(dexseqResults))[1]]
  }else{
    y <- dexseqResults[,foldChangeColumn]
  }
  df <- data.frame( x, y, object$padj < alpha )
#  ylim=c(-2, 2)
  plotMA(df, ylim, ... )
}

setMethod("plotMA", signature(object="DEXSeqResults"),
          plotMA.DEXSeqResults)

show.DEXSeqResults <- function(object){
  cat("\n")
  cat(mcols(object)[colnames(object) == "pvalue","description"])
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


subsetByOverlaps.DEXSeqResults <- function( x, ranges, maxgap = -1L, minoverlap = 0L,
         type = c("any", "start", "end", "within", "equal"),
         ignore.strand = FALSE ){
  stopifnot( is( x, "DEXSeqResults") )
  genomicData <- x$genomicData
  overlaps <- findOverlaps( query=genomicData, subject=ranges, maxgap=maxgap,
    minoverlap=minoverlap, type=type,
    ignore.strand=ignore.strand )
  x[queryHits( overlaps ),]
}

setMethod("subsetByOverlaps", signature(x="DEXSeqResults", ranges="GenomicRanges"),
          subsetByOverlaps.DEXSeqResults)

findOverlaps.DEXSeqResults <- function( query, subject, maxgap = -1L, minoverlap = 0L,
         type = c("any", "start", "end", "within", "equal"),
         ignore.strand = FALSE ){
  stopifnot( is( query, "DEXSeqResults") )
  genomicData <- query$genomicData
  overlaps <- findOverlaps( query=genomicData, subject=subject, maxgap=maxgap,
    minoverlap=minoverlap, type=type,
    ignore.strand=ignore.strand)
  overlaps
}

setMethod("findOverlaps", signature(query="DEXSeqResults", subject="GenomicRanges"),
          findOverlaps.DEXSeqResults)

setMethod("[", "DEXSeqDataSet", function(x, i, j, ..., drop = FALSE) {
    x <- callNextMethod()
    colData(x) <- droplevels(colData(x))
    x@modelFrameBM <- makeBigModelFrame(x)
    validObject(x)
    x
})

setReplaceMethod("colData", c("DEXSeqDataSet", "DataFrame"),
    function(x, ..., value)
{
    x <- callNextMethod()
    x@modelFrameBM <- makeBigModelFrame(x)
    validObject(x)
    x
})

setReplaceMethod("$", "DEXSeqDataSet",
    function(x, name, value)
{
    x <- callNextMethod()
    x@modelFrameBM <- makeBigModelFrame(x)
    validObject(x)
    x
})
