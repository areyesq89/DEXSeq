
prepareAnnotationForDEXSeq <- function( transcriptDb, aggregateGenes=FALSE, includeTranscripts=TRUE )
{

  stopifnot( is( transcriptDb, "TranscriptDb" ) )
  exonsByGene <- exonsBy(transcriptDb, by="gene")
  
  exonicParts <- disjoin( unlist(exonsByGene) )

  if( !aggregateGenes ){
     overlaps <- findOverlaps( exonicParts, exonsByGene )
     geneNames <- names(exonsByGene)[ subjectHits( overlaps ) ]
     aggregateGeneNames <- split( geneNames, queryHits(overlaps) )
     toRemove <- names(aggregateGeneNames)[which( sapply( aggregateGeneNames, length ) > 1 )]
     toRemove <- as.numeric( toRemove )
     if( length( toRemove ) > 0 ){
       exonicParts <- exonicParts[-toRemove]
       geneNames <- aggregateGeneNames[-toRemove]
     }
     mcols( exonicParts )$geneNames <- unlist( geneNames )
  }else{ 
     foGG <- findOverlaps(exonsByGene, exonsByGene)
     splitByGene <- split(subjectHits(foGG), queryHits(foGG))
     aggregateGeneNames <- sapply( splitByGene, function(i){
                             paste(names(exonsByGene)[i],collapse="+")} )
     foEG <- findOverlaps(exonicParts, exonsByGene, select="first")
     mcols(exonicParts)$geneNames <- aggregateGeneNames[foEG]
  }

  if( includeTranscripts ){
     exonsByTranscript <- exonsBy( transcriptDb, by="tx", use.names=TRUE )
     foET <- findOverlaps(exonicParts, exonsByTranscript)
     splitByExonicPart <- split(subjectHits(foET), queryHits(foET))
     mcols(exonicParts)$transcripts <- sapply(splitByExonicPart, function(i){
        paste(names(exonsByTranscript)[i],collapse=";")})
   }

   exonicParts <- exonicParts[order(mcols(exonicParts)$geneNames)]
   mcols(exonicParts)$exonic_part_number <- do.call( c, lapply(
      split(mcols(exonicParts)$geneNames,mcols(exonicParts)$geneNames), 
      function(z){ seq(along=z) } ) )


   mcols( exonicParts )$exonID <- sprintf("E%03.0f", mcols( exonicParts )$exonic_part_number )
   exonicParts

}

countReadsForDEXSeq <- function( exonicParts, bamFileList, scanBamParam=ScanBamParam(), singleEnd=TRUE, ignoreStrand=TRUE, 
   mode=function(reads, features, ignore.strand){ countOverlaps( features, reads, ignore.strand=ignoreStrand)} )
{
   stopifnot( is( bamFileList, "BamFileList" ) )
   stopifnot( is( exonicParts, "GRanges" ) )

   exonHits <- summarizeOverlaps( exonicParts, bamFileList, mode=mode, singleEnd=singleEnd, ignore.strand=ignoreStrand, param=scanBamParam)
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



