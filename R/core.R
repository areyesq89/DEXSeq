rmDepCols <- function(m)
{
    q <- qr(m)
    if (q$rank < ncol(m)) 
        m[, -q$pivot[(q$rank + 1):ncol(m)]]
    else m
}


#BPPARAM <- MulticoreParam(workers=1)
#fullModel <- design(object)
#reducedModel <- ~ sample + exon
vst <- function(x,  object ){
  if( is( object, "DEXSeqDataSet") ){
     coefs <- attr(dispersionFunction(object), "coefficients")
  }
  if( is( object, "DEXSeqResults" ) ){
     coefs <- object@coefs
  }
    (2/(sqrt(coefs["asymptDisp"]))) * log(2 * coefs["asymptDisp"] * 
      sqrt(x) + 2 * sqrt(coefs["asymptDisp"] * (coefs["extraPois"] + 
        1 + coefs["asymptDisp"] * x))) - (2/(sqrt(coefs["asymptDisp"]))) * 
          log(2 * sqrt(coefs["asymptDisp"] * (coefs["extraPois"] + 1)))
}


testForDEU <-
  function(object,
           fullModel=design(object),
           reducedModel= ~ sample + exon,
           BPPARAM=MulticoreParam(workers=1))
{
  if (is.null(sizeFactors(object)) & is.null(normalizationFactors(object))) {
    stop("first call estimateSizeFactors or provide a normalizationFactor matrix")
  }
  if (is.null(dispersions(object))) {
    stop("first call estimateDispersions")
  }

  splitParts <- sort(
    rep(seq_len(BPPARAM$workers), 
    length.out=nrow(object) ) )
  splitObject <- split( object, splitParts )

  splitObject <- bplapply( splitObject,
    function(x){
      x <- nbinomLRT( x, reduced = reducedModel, full=fullModel )
    }, BPPARAM=BPPARAM )

  mcols( object ) <- mcols( do.call(rbind, splitObject) )
  assays( object ) <- assays( do.call(rbind, splitObject) )

  extraAttributes <- setdiff( names( attributes(splitObject[[1]]) ),  names( attributes(object) ) )

  for( atr in extraAttributes ){
    attr( object, atr ) <- attr( splitObject[[1]], atr )
  }

  object
 
}

estimateExonFoldChanges <- function( object,
                                    fitExpToVar = "condition",
                                    denominator = "",
                                    BPPARAM=MulticoreParam(workers=1) )
{
    stopifnot(is(object, "DEXSeqDataSet"))
    if (any(is.na(sizeFactors(object)))) {
      stop("Please estimate sizeFactors first\n")
    }
    if (!fitExpToVar %in% colnames(sampleAnnotation(object))) {
      stop(sprintf("%s parameter is not in the colData", fitExpToVar))
    }
    if ( is.null( dispersions(object) ) ){
      stop("please call estimateDispersions first")
    }
    if ( length( attr( object, "test" ) ) == 0 ){
      stop("please call testForDEU first")
    }
    frm <- as.formula(paste("count ~", fitExpToVar, "* exon"))
    testablegenes <- as.character( unique( groupIDs(object)[!is.na( results(object)$padj )] ) )
    groups <- groupIDs(object)
    features <- featureIDs(object)
    exonCounts <- featureCounts( object )
    disps <- dispersions(object)
    maxMf <- object@modelFrameBM
    rowsPerSample <- split(seq_len(nrow(maxMf)), maxMf$sample)
    geteffects <- function(geneID) {
#        print( geneID )
        rows <- groups %in% geneID
        numExons <- sum(rows)
        newMf <- maxMf[as.vector( sapply(rowsPerSample, "[", seq_len(numExons)) ),]
        newMf <- droplevels( newMf )
        newMf$count <- as.vector( exonCounts[rows,] )
        newMf$dispersion <- rep( disps[rows], ncol(exonCounts) )
        coefficients <- fitAndArrangeCoefs( frm, balanceExons = TRUE, mf=newMf)
        if (is.null(coefficients)) {
            return(coefficients)
        }
        ret <- t(getEffectsForPlotting(coefficients, averageOutExpression = TRUE, 
            groupingVar = fitExpToVar))
        rownames(ret) <- paste(geneID, rownames(ret), sep = ":")
        return(ret)
    }
   
    alleffects <- bplapply( testablegenes, geteffects, BPPARAM=BPPARAM )
    alleffects <- do.call(rbind, alleffects)
    alleffects <- vst(exp(alleffects), object)
    toadd <- matrix(NA, nrow = nrow(object), ncol = ncol(alleffects))
    rownames(toadd) <- rownames(object)
    colnames(toadd) <- colnames(alleffects)
    toadd[rownames(alleffects), colnames(alleffects)] <- alleffects
    toadd <- DataFrame(toadd)
    elementMetadata(toadd) <- DataFrame(
                                        type=rep("DEXSeq results",
                                          ncol(toadd)),
                                        description=rep("exon usage coefficient",
                                          ncol(toadd) ) )
    toadd2 <- matrix(NA, nrow = nrow(object), ncol = ncol(alleffects) )
    if (denominator == "") {
        denominator <- as.character(sampleAnnotation(object)[[fitExpToVar]][1])
    }
    colRemove <- colnames(alleffects) %in% denominator
    stopifnot(any(colRemove))
    denoCol <- which(colnames(alleffects) == denominator)
    alleffects <- log2(alleffects/alleffects[, denoCol])
    colnames(alleffects) <- sprintf("log2fold_%s_%s", colnames(alleffects), 
            denominator)
    colnames(toadd2) <- colnames(alleffects)
    rownames(toadd2) <- rownames(object)
    toadd2[rownames(alleffects), colnames(alleffects)] <- alleffects
    toadd2 <- toadd2[,-denoCol, drop=FALSE]
    toadd2 <- DataFrame(toadd2)
    elementMetadata(toadd2) <- DataFrame(type=rep("DEXSeq results",
                                          ncol(toadd2)),
                                       description=rep("relative exon usage fold change",
                                          ncol(toadd2) ) )
    allAdd <- cbind(toadd, toadd2)
    mcols(object) <- cbind(mcols(object), allAdd)
    object
}


DEXSeqResults <- function( object ){
  stopifnot( is(object, "DEXSeqDataSet"))
  LRTresults <- results(object)
  LRTresults$exonBaseMean <- rowMeans(featureCounts(object))
  LRTresults$featureID <- mcols(object)$featureID
  LRTresults$groupID <- mcols(object)$groupID
  LRTresults$dispersion <- mcols(object)$dispersion
  
  LRTresults <- LRTresults[,c("groupID", "featureID", "exonBaseMean", "dispersion", "stat", "pvalue", "padj")]
  elementMetadata( LRTresults )[colnames(LRTresults) %in% c("groupID", "featureID", "exonBaseMean"),"type"] <- "input"
  elementMetadata( LRTresults )[colnames(LRTresults) %in% "groupID","description"] <- "group/gene identifier"
  elementMetadata( LRTresults )[colnames(LRTresults) %in% "featureID","description"] <- "feature/exon identifier"
  elementMetadata( LRTresults )[colnames(LRTresults) %in% "exonBaseMean","description"] <- "mean of the counts across samples in each feature/exon"
  elementMetadata( LRTresults )[colnames(LRTresults) %in% "dispersion","description"] <- "exon dispersion estimate"
  toadd <- mcols(object)[,elementMetadata( mcols(object ) )$type == "DEXSeq results", drop=FALSE]
  LRTresults <- cbind( LRTresults, toadd )
  genomicData <- rowData(object)
  mcols(genomicData) <- NULL
  LRTresults$genomicData <- genomicData
  LRTresults$countData <- featureCounts(object)
  LRTresults$transcripts <- mcols(object)$transcripts
  elementMetadata( LRTresults )[colnames(LRTresults) %in% c("genomicData", "countData", "transcripts"),"type"] <- "input"
  elementMetadata( LRTresults )[colnames(LRTresults) %in% "genomicData","description"] <- "GRanges object of the coordinates of the exon/feature"
  elementMetadata( LRTresults )[colnames(LRTresults) %in% "countData","description"] <- "matrix of integer counts, of each column containing a sample"
  elementMetadata( LRTresults )[colnames(LRTresults) %in% "transcripts","description"] <- "list of transcripts overlapping with the exon"
  dxr <-
      new("DEXSeqResults",
          LRTresults,
          modelFrameBM=object@modelFrameBM,
          sampleData = sampleAnnotation(object),
          coefs = attr(dispersionFunction(object), "coefficients") )

  dxr
}

DEXSeq <- function( object,
                   fullModel=design(object),
                   reducedModel = ~ sample + exon,
                   BPPARAM=MulticoreParam(workers=1), fitExpToVar="condition", quiet=TRUE ){
  stopifnot(is( object, "DEXSeqDataSet") )
  object <- estimateSizeFactors( object )
  object <- estimateDispersions( object, formula=fullModel, BPPARAM=BPPARAM, quiet=TRUE)
  object <- testForDEU( object, reducedModel=reducedModel, fullModel=fullModel, BPPARAM=BPPARAM )
  object <- estimateExonFoldChanges( object, fitExpToVar=fitExpToVar )
  res <- DEXSeqResults( object )
  res
}
