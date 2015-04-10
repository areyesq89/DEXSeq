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
vst <- function( x,  object ){
  if( attr(object@dispersionFunction, "fitType") != "parametric" ){
    warnings("The dispersion function is not parametric, the DEXSeq vst won't be applied to the data\n")
    return( x )
  }
  coefs <- attr(object@dispersionFunction, "coefficients")

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

  allVars <- all.vars(reducedModel)
  if( any(!allVars %in% colnames( colData(object) )) ){
     notPresent <- allVars[!allVars %in% colnames( colData(object) )]
     notPresent <- paste(notPresent, collapse=",")
     stop(sprintf("the variables '%s' of the parameter 'reducedModel' are not specified in the columns of colData", notPresent ) )
  }

  allVars <- all.vars(fullModel)
  if( any(!allVars %in% colnames( colData(object) )) ){
     notPresent <- allVars[!allVars %in% colnames( colData(object) )]
     notPresent <- paste(notPresent, collapse=",")
     stop(sprintf("the variables '%s' of the parameter 'fullModel' are not specified in the columns of colData", notPresent ) )
  }

  fullModelMatrix <- 
    rmDepCols( model.matrix( fullModel, as.data.frame(colData(object)) ) )

  reducedModelMatrix <- 
    rmDepCols( model.matrix( reducedModel, as.data.frame(colData(object)) ) )

  splitParts <- sort(
    rep(seq_len(BPPARAM$workers), 
    length.out=nrow(object) ) )
  splitObject <- split( object, splitParts )

  splitObject <- bplapply( splitObject,
    function(x){
      x <- nbinomLRT( x, reduced = reducedModelMatrix, full=fullModelMatrix )
    }, BPPARAM=BPPARAM )

  mergeObject <- do.call(rbind, splitObject)
  matchedNames <- match( rownames(object), rownames(mergeObject))  
  mcols(object) <- mcols( mergeObject )[matchedNames,]
  assays(object) <- assays(mergeObject[matchedNames,])

  extraAttributes <- setdiff( names( attributes(splitObject[[1]]) ),  names( attributes(object) ) )

  for( atr in extraAttributes ){
    attr( object, atr ) <- attr( splitObject[[1]], atr )
  }

  object
 
}

estimateExonFoldChanges <- function( object,
                                    fitExpToVar = "condition",
                                    denominator = "",
                                    BPPARAM=MulticoreParam(workers=1), 
                                    maxRowsMF=3000)
{
    stopifnot(is(object, "DEXSeqDataSet"))
    if (any(is.na(sizeFactors(object)))) {
      stop("Please estimate sizeFactors first\n")
    }
    if (!fitExpToVar %in% colnames(sampleAnnotation(object))) {
      stop(sprintf("The value of the parameter fitExpToVar,'%s', is not a column name of colData", fitExpToVar))
    }
    if ( is.null( dispersions(object) ) ){
      stop("please call estimateDispersions first")
    }
    if ( length( attr( object, "test" ) ) == 0 ){
      stop("please call testForDEU first")
    }
    frm <- as.formula(paste("count ~", fitExpToVar, "* exon"))
    notNAs <- !is.na( results(object, 
       filter=rowMeans( featureCounts(object, normalized=TRUE) ))$padj )
    testablegenes <- unique(groupIDs(object)[notNAs])
    groups <- groupIDs(object)
    features <- featureIDs(object)
    exonCounts <- featureCounts( object )
    disps <- dispersions(object)
    disps[is.na( disps )] <- 1e-8
    mf <- object@modelFrameBM
    rowsPerSample <- split(seq_len(nrow(mf)), mf$sample)
    numsamples <- nrow( sampleAnnotation(object) )
    features <- featureIDs(object)
    countsAll <- featureCounts(object)
    geteffects <- function(geneID){
#        print( geneID )
        rt <- groups %in% geneID
        numexons <- sum(rt)
        newMf <- mf[as.vector( sapply( split( seq_len(nrow(mf)), mf$sample ), "[", seq_len( numexons ) ) ),]
        featuresInGene <- features[rt]
        newMf$exon <- factor( rep( featuresInGene, numsamples ) )
        countsThis <- countsAll[rt,]
        rownames(countsThis) <- gsub("\\S+:", "", rownames(countsThis))
        dispsThis <- disps[rt]
        names(dispsThis) <- features[rt]
        for( i in seq_len(nrow(newMf))){
           newMf[i,"dispersion"] <- dispsThis[as.character(newMf[i,"exon"])]
           newMf[i,"count"] <- countsThis[as.character(newMf[i,"exon"]), as.character(newMf[i,"sample"])]
        }
        newMf <- droplevels( newMf )
        coefficients <- fitAndArrangeCoefs( frm, balanceExons = TRUE, mf=newMf, maxRowsMF=maxRowsMF)
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
  LRTresults <- results(object, filter=rowMeans( featureCounts(object, normalized=TRUE) ) )
  LRTresults$exonBaseMean <- rowMeans(featureCounts(object, normalized=TRUE))
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
  genomicData <- rowRanges(object)
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
          dispersionFunction = object@dispersionFunction )

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
