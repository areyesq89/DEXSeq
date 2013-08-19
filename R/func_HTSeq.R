read.HTSeqCounts <- function( countfiles, design, flattenedfile=NULL )
{
   lf <- lapply( countfiles, function(x)
      read.table( x, header=FALSE,stringsAsFactors=FALSE ) )
   if( !all( sapply( lf[-1], function(x) all( x$V1 == lf[1]$V1 ) ) ) )
      stop( "Count files have differing gene ID column." )
   dcounts <- sapply( lf, `[[`, "V2" )
   rownames(dcounts) <- lf[[1]][,1]
   rownames(dcounts)
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
      transcripts <- gsub("\\+", ";", transcripts)
      exonids <- gsub(".*exonic_part_number\\s(\\S+).*", "\\1", aggregates$attr)
      exoninfo<-data.frame(chr=aggregates$chr, start=aggregates$start, end=aggregates$end, strand=aggregates$strand)
      rownames( exoninfo ) <- paste( aggregates$gene_id, exonids, sep=":E" )
      names(transcripts) <- rownames(exoninfo)
      if (!all( rownames(dcounts) %in% rownames(exoninfo) )){
         stop("Count files do not correspond to the flattened annotation file")
      }
      matching <- match(rownames(dcounts), rownames(exoninfo))
      ecs<-newExonCountSet(countData=dcounts, design=design, geneIDs=genesrle, exonIDs=exons, exonIntervals=exoninfo[matching,], transcripts=transcripts[matching])
      ecs@annotationFile <- flattenedfile
      pData(ecs)$countfiles <- countfiles
      return(ecs)
   }else{
      return(newExonCountSet(countData=dcounts, design=design, geneIDs=genesrle, exonIDs=exons))
   }
}
