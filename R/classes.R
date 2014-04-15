setClass("DEXSeqDataSet",
   contains = "DESeqDataSet",
         representation = representation( modelFrameBM = "data.frame") )

DEXSeqDataSet <- function( countData, sampleData, design= ~ sample + exon + condition:exon , featureID, groupID, featureRanges=NULL, transcripts=NULL)
{

  stopifnot( class( countData ) %in% c("matrix", "data.frame"))
  countData <- as.matrix( countData )
  stopifnot( class( featureID ) %in% c("character", "factor"))
  stopifnot( class( groupID ) %in% c("character", "factor"))
  stopifnot( class( sampleData ) %in% c("data.frame"))

  modelFrame <- cbind(
    sample = rownames(sampleData), sampleData )
  modelFrame <- rbind( cbind(modelFrame, exon = "this"), 
    cbind(modelFrame, exon = "others"))
  rownames(modelFrame) <- NULL
  colData <- DataFrame( modelFrame )

  if( !"exon" %in% all.vars( design ) ){
    stop("The formula does not specify a contrast with the variable 'exon'")
  }
  
  forCycle <- split( 1:nrow( countData ), as.character( groupID ) )
  others <- lapply( forCycle, function(i){
    sct <- countData[i, , drop = FALSE]
    rs <- t(sapply(1:nrow(sct), function(r) colSums(sct[-r, , drop = FALSE])))
    rownames(rs) <- rownames(sct)
    rs
  })

#  others <- tapply(1:nrow(countData), as.character(groupID), function(i) {
#    sct <- countData[i, , drop = FALSE]
#    t(sapply(1:nrow(sct), function(r) colSums(sct[-r, , drop = FALSE])))
#  })

  others <- do.call(rbind, others)
  stopifnot( all( rownames(countData) %in% rownames(others) ) )
  others <- others[rownames(countData),]
  nCountData <- cbind( countData, others ) 

  if( !is.null(featureRanges) ){ 
    stopifnot( class(featureRanges) %in% c("GRanges", "GRangesList"))
    se <- SummarizedExperiment( nCountData, colData=colData, rowData=featureRanges )
  }else{
    se <- SummarizedExperiment( nCountData, colData=colData )
  }

  names(assays(se))[1] = "counts"
  mcols( se )$featureID <- featureID
  mcols( se )$groupID <- groupID
  mcols( se )$exonBaseMean <- rowMeans( countData )
  mcols( se )$exonBaseVar <- rowVars( countData )
  
  if( !is.null(transcripts) ){
    mcols(se)$transcripts <- transcripts
  }
     
  rownames(se) <- paste( groupID, featureID, sep=":")
  
  dds <- DESeqDataSet( se, design, ignoreRank=TRUE )
  
  maxGene <- names(which.max(table(groupID)))
  rows <- mcols(dds)$groupID %in%  maxGene
  numExons <- sum( rows )
  
  exonCol <-
    rep(factor(featureID[rows]), nrow(sampleData))

  modelFrame <- data.frame(
    sample=rep( rownames(sampleData), each=numExons),
    exon = exonCol )
  
  varNames <- colnames( sampleData )
  for( i in varNames ){
    modelFrame[[i]] <- rep( sampleData[[i]], each=numExons )
  }

  modelFrame$dispersion <- NA
  modelFrame$sizeFactor <- NA
  modelFrame$count <- NA
  
  dxd <- new( "DEXSeqDataSet", dds, modelFrameBM=modelFrame )

  return(dxd)
}


setValidity( "DEXSeqDataSet", function( object ) {
  stopifnot(
      c("sample", "exon", "dispersion", "sizeFactor", "count")
            %in% colnames( object@modelFrameBM ) )
  TRUE
} )

setClass("DEXSeqResults",
   contains = "DataFrame",
         representation = representation( modelFrameBM = "data.frame", sampleData="DataFrame", dispersionFunction = "function") )


setValidity( "DEXSeqResults", function( object ){
    stopifnot( "sample" %in% colnames( object@sampleData ) )
    stopifnot( colnames(object$countData) == as.character(object@sampleData$sample) )    
    TRUE
})

###########################
#### ACCESSOR FUNCTIONS####
###########################


featureCounts <- function( object, normalized=FALSE ){
  validObject(object)
  res <- counts(object, normalized=normalized)[,colData(object)$exon == "this"]
  colnames( res ) <- sampleAnnotation(object)$sample
  res
}

featureIDs <- function(object){
  validObject(object)
  mcols(object)$featureID
}

`featureIDs<-` <- function( object, value ) {
   stopifnot( is( object, "DEXSeqDataSet" ) )
   mcols(object)$featureID <- value
   rownames(object) <- paste( mcols(object)$groupID, mcols(object)$featureID, sep=":" )
   validObject(object)
   object
}

exonIDs <- function(object){
  validObject(object)
  featureIDs(object)
}

`exonIDs<-` <- function( object, value ) {
  object <- `featureIDs<-`( object, value )
  object
}

groupIDs <- function( object ){
  validObject( object )
  mcols( object )$groupID
}

`groupIDs<-` <- function( object, value ) {
   stopifnot( is( object, "DEXSeqDataSet" ) )
   mcols( object )$groupID <- value
   rownames(object) <- paste( mcols(object)$groupID, mcols(object)$featureID, sep=":" )
   validObject( object )
   object
}

geneIDs <- function( object ){
  validObject( object )
  groupIDs( object )
}

`geneIDs<-` <- function( object, value ) {
  object <- `groupIDs<-`(object, value)
  object
}


sampleAnnotation <- function( object ){
  validObject( object )
  colData( object )[colData( object )$exon == "this",!colnames(colData( object )) %in% "exon"]
}

#################
###FROM HTSEQ####
#################


DEXSeqDataSetFromHTSeq <- function( countfiles, sampleData, design= ~ sample + exon + condition:exon, flattenedfile=NULL )
{
   if( !all( sapply(countfiles, class) == 'character' ) ){
      stop("The countfiles parameter must be a character vector")
   }
   lf <- lapply( countfiles, function(x)
      read.table( x, header=FALSE,stringsAsFactors=FALSE ) )
   if( !all( sapply( lf[-1], function(x) all( x$V1 == lf[1]$V1 ) ) ) )
      stop( "Count files have differing gene ID column." )
   dcounts <- sapply( lf, `[[`, "V2" )
   rownames(dcounts) <- lf[[1]][,1]
   dcounts <- dcounts[ substr(rownames(dcounts),1,1)!="_", ]
   rownames(dcounts) <- sub(":", ":E", rownames(dcounts))
   colnames(dcounts) <- countfiles
   splitted <- strsplit(rownames(dcounts), ":")
   exons <- sapply(splitted, "[[", 2)
   genesrle <- sapply( splitted, "[[", 1)
   if(!is.null(flattenedfile)){
      aggregates<-read.delim(flattenedfile, stringsAsFactors=FALSE, header=FALSE)
      colnames(aggregates)<-c("chr", "source", "class", "start", "end", "ex", "strand", "ex2", "attr")
      aggregates<-aggregates[which(aggregates$class =="exonic_part"),]
      aggregates$attr <- gsub("\"|=|;", "", aggregates$attr)
      aggregates$gene_id <- sub(".*gene_id\\s(\\S+).*", "\\1", aggregates$attr)
      transcripts <- gsub(".*transcripts\\s(\\S+).*", "\\1", aggregates$attr)
      transcripts <- strsplit(transcripts, "\\+")
      exonids <- gsub(".*exonic_part_number\\s(\\S+).*", "\\1", aggregates$attr)
      exoninfo<-GRanges(
                        as.character(aggregates$chr),
                        IRanges(start=aggregates$start, end=aggregates$end),
                        strand=aggregates$strand)
      names( exoninfo ) <- paste( aggregates$gene_id, exonids, sep=":E" )
      names(transcripts) <- rownames(exoninfo)
      if (!all( rownames(dcounts) %in% names(exoninfo) )){
         stop("Count files do not correspond to the flattened annotation file")
      }
      matching <- match(rownames(dcounts), names(exoninfo))
      stopifnot( all( names( exoninfo[matching] ) == rownames(dcounts) ) )
      stopifnot( all( names( transcripts[matching] ) == rownames(dcounts) ) )
      
      dxd <- DEXSeqDataSet( dcounts, sampleData, design, exons, genesrle,
                           exoninfo[matching], transcripts[matching] )
      return(dxd)
   }else{
      dxd <- DEXSeqDataSet( dcounts, sampleData, design, exons, genesrle)
      return(dxd)
   }
}

DEXSeqDataSetFromSE <- function( SE, design= ~ sample + exon + condition:exon ){
  if( !all( c("gene_id", "tx_name", "exonic_part") %in% colnames( mcols( SE  ) ) ) ){
    stop("make sure your SummarizedExperiment object contain the columns gene_id, tx_name and exonic_part")
  }
  groupID <- as.character( mcols(SE)$gene_id )
  featureID <- sprintf("E%3.3d",  mcols(SE)$exonic_part )
  transcripts <- as.list( mcols(SE)$tx_name )
  sampleData <- as.data.frame( colData(SE) )
  design <- design
  mcols(SE) <- NULL
  featureRanges <- rowData(SE)
  countData <- assay(SE)
  dxd <- DEXSeqDataSet( countData,
      sampleData,
      design,
      featureID,
      groupID,
      featureRanges=featureRanges,
      transcripts=transcripts )
  dxd
}
