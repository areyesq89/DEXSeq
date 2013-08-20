prepareAnnotationForDEXSeq <- function(transcriptDb, aggregateGenes = FALSE,
                 includeTranscripts = TRUE)
{
    .Deprecated("disjointExons", package="DEXSeq")
}

countReadsForDEXSeq <- function (exonicParts, bamFileList, scanBamParam = ScanBamParam(), 
        singleEnd = TRUE, ignoreStrand = TRUE, mode = function(features, reads, 
        ignore.strand, inter.feature=FALSE) {
        countOverlaps(features, reads, ignore.strand = ignoreStrand)
    })
{
    .Deprecated("summarizedOverlaps", package="DEXSeq")
}



buildExonCountSet <- function( summarizedExperiment, design, exonicParts ){
   stopifnot( is( exonicParts, "GRanges" ) )
   stopifnot( is( summarizedExperiment, "SummarizedExperiment" ) )
   stopifnot( all( colnames( mcols( exonicParts ) ) %in% c("gene_id", "tx_name", "exonic_part") ))
   chr <- rep( as.character(seqnames(exonicParts)@values), times=seqnames(exonicParts)@lengths)
   strand <- rep( as.character( strand(exonicParts)@values), times=strand(exonicParts)@lengths)
   geneIDsReady <- sapply( as.list(mcols(exonicParts)$gene_id), paste, collapse="+" )
   transcriptsReady <- sapply( as.list(mcols(exonicParts)$tx_name), paste, collapse=";" )
   newExonCountSet( 
      assay( summarizedExperiment), 
      design=design, 
      geneIDs=geneIDsReady, 
      exonIDs=sprintf("E%3.3d", unlist( mcols( exonicParts )$exonic_part)),
      transcripts=transcriptsReady,
      exonIntervals= data.frame(chr=chr, start=start(exonicParts), end=end(exonicParts), strand=strand)
   )
}



