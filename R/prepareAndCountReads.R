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
    stopifnot(is(bamFileList, "BamFileList"))
    stopifnot(is(exonicParts, "GRanges"))
    exonHits <- summarizeOverlaps(exonicParts, bamFileList, mode = mode, 
        singleEnd = singleEnd, fragments=FALSE, ignore.strand = ignoreStrand, 
        param = scanBamParam, inter.feature=FALSE)
    exonHits
}



buildExonCountSet <- function( summarizedExperiment, design, exonicParts ){
   stopifnot( is( exonicParts, "GRanges" ) )
   stopifnot( is( summarizedExperiment, "SummarizedExperiment" ) )
   stopifnot( length( mcols(exonicParts)$geneNames ) == length( exonicParts ) )
   stopifnot( length( mcols(exonicParts)$exonID ) == length( exonicParts ) )
   chr <- rep( seqnames(exonicParts)@values, times=seqnames(exonicParts)@lengths)
   strand <- rep( strand(exonicParts)@values, times=strand(exonicParts)@lengths)
   newExonCountSet( 
      assay( summarizedExperiment), 
      design=design, 
      geneIDs=mcols( exonicParts )$geneNames, 
      exonIDs=mcols( exonicParts )$exonID,
      transcripts=mcols( exonicParts )$transcripts,
      exonIntervals= data.frame(chr=chr, start=start(exonicParts), end=end(exonicParts), strand=strand)
   )
}



