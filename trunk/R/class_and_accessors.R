setClass( "ExonCountSet", 
   contains = "eSet",
   representation = representation( 
      designColumns = "character",
      dispFitCoefs = "numeric"
   ),
   prototype = prototype( new( "VersionedBiobase",
      versions = c( classVersion("eSet"), ExonCountSet = "1.0.1" ) ) )
)


newExonCountSet <- function( countData, design, geneIDs, exonIDs, exonIntervals=NULL, transcripts=NULL){

   countData <- as.matrix( countData )
   if( any( round( countData ) != countData ) ){
      stop( "The countData is not integer." )}
   mode( countData ) <- "integer"

   if( is( design, "matrix" ) ){
      design <- as.data.frame( design )}

#   if(!(is(design, "data.frame") || is(design, "AnnotatedDataFrame"))){
 #     countData <- countData[,order(as.character(design))]
  #    design <- design[order(as.character(design))]
 #  }else{
 #     countData <- countData[,order(design$condition)]
 #     design <- design[order(design$condition),]
 #  }

   phenoData <- annotatedDataFrameFrom( countData, byrow=FALSE )
   featureData <- annotatedDataFrameFrom( countData, byrow=TRUE )
      
   phenoData$sizeFactor <- rep( NA_real_, ncol(countData) )
   varMetadata( phenoData )[ "sizeFactor", "labelDescription" ] <- "size factor (relative estimate of sequencing depth)"

   geneIDs <- as.factor( geneIDs )
   if( length(geneIDs) != nrow(countData) )
      stop( "geneIDs must be of the same length as the number of columns in countData")

   featureData$geneID <- geneIDs
   varMetadata( featureData )[ "geneID", "labelDescription" ] <- "ID of gene to which the exon belongs"

   exonIDs <- as.character( exonIDs )
   if( length(exonIDs) != nrow(countData) )
      stop( "exonIDs must be of the same length as the number of columns in countData")

   featureData$exonID <- exonIDs
   varMetadata( featureData )[ "exonID", "labelDescription" ] <- "exon ID (unique only within a gene)"

   if( is.null(exonIntervals) ){
      exonIntervals <- data.frame(
         chr    = rep( NA_character_, nrow( countData ) ), 
         start  = rep( NA_integer_,   nrow( countData ) ), 
         end    = rep( NA_integer_,   nrow( countData ) ), 
         strand = rep( NA_character_, nrow( countData ) ) ) }

   featureData$dispersion_CR_est <- rep( NA_real_, nrow( countData ) )
   varMetadata( featureData )[ "dispersion_CR_est", "labelDescription" ] <- "exon dispersion (Cox-Reid estimate)"

   featureData$dispersion <- rep( NA_real_, nrow( countData ) )
   varMetadata( featureData )[ "dispersion", "labelDescription" ] <- "exon dispersion (value used in test)"

   featureData$pvalue <- rep( NA_real_, nrow( countData ) )
   varMetadata( featureData )[ "pvalue", "labelDescription" ] <- "p-value from testForDEU"

   featureData$padjust <- rep( NA_real_, nrow( countData ) )
   varMetadata( featureData )[ "padjust", "labelDescription" ] <- "BH adjusted p-value"

   exonIntervals <- as.data.frame( exonIntervals )
   
   # in case it was a GRanges object before, change the colname:
   if( "seqnames" %in% colnames(exonIntervals) ){
      colnames(exonIntervals)[ colnames(exonIntervals) == "seqnames" ] <- "chr"   }
   
   if( !all( c( "chr", "start", "end", "strand" ) %in% colnames(exonIntervals) ) ){
      stop( "exonIntervals must be a data frame with columns 'chr', 'start', 'end', and 'strand'." )}

   if(is.null(transcripts)){
      transcripts <- rep(NA_character_, nrow( countData ) )}

   if(!is(transcripts, "character")){
      stop("transcript information must be a character vector")}
   
   featureData$chr    <- factor( exonIntervals$chr )
   featureData$start  <- exonIntervals$start
   featureData$end    <- exonIntervals$end
   featureData$strand <- factor( exonIntervals$strand )
   featureData$transcripts <- transcripts
   varMetadata( featureData )[ "chr",    "labelDescription" ] <- "chromosome of exon"
   varMetadata( featureData )[ "start",  "labelDescription" ] <- "start of exon"
   varMetadata( featureData )[ "end",    "labelDescription" ] <- "end of exon"
   varMetadata( featureData )[ "strand", "labelDescription" ] <- "strand of exon"
   varMetadata( featureData )[ "transcripts", "labelDescription" ] <- "transcripts in which this exon is contained"
   
   if( is( design, "data.frame" ) || is( design, "AnnotatedDataFrame" ) ) {
      stopifnot( nrow( design ) == ncol( countData ) )
      design <- as( design, "AnnotatedDataFrame" )
      dimLabels(design) <- dimLabels(phenoData)
      rownames( pData(design) ) <- rownames( pData(phenoData) )
      phenoData <- combine( phenoData, design )
      rvft <- c( `_all` = NA_character_ )
      designColumns <- varLabels(design)
   } else {
      design <- factor( design, levels=unique(design))
      stopifnot( length( design ) == ncol( countData ) )
      phenoData$`condition` <- factor( design )
      varMetadata( phenoData )[ "condition", "labelDescription" ] <- "experimental condition, treatment or phenotype"
      designColumns <- "condition"
   }
   ecs <- new( "ExonCountSet",
      assayData = assayDataNew( "environment", counts=countData ),
      phenoData = phenoData, 
      featureData = featureData,
      designColumns = designColumns,
      dispFitCoefs = c( NA_real_, NA_real_ )
      )
   ecs
}

setValidity( "ExonCountSet", function( object ) {

   if( !all( object@designColumns %in% names(pData(object)) ) )
      return( "Not all designColumns appear in phenoData." )

   if( ! "sizeFactor" %in% names(pData(object)) )
      return( "phenoData does not contain a 'sizeFactor' column.")
   if( ! is( pData(object)$`sizeFactor`, "numeric" ) )
      return( "The 'sizeFactor' column in phenoData is not numeric." )

   if( ! "geneID" %in% names(fData(object)) )
      return( "featureData does not contain a 'geneID' column.")
   if( ! is( fData(object)$geneID, "factor" ) )
      return( "The 'geneID' column in fData is not a factor." )

   if( ! "exonID" %in% names(fData(object)) )
      return( "featureData does not contain an 'exonID' column.")
   if( ! is( fData(object)$exonID, "character" ) )
      return( "The 'exonID' column in fData is not a character vector." )

   if( ! "chr"  %in% names(fData(object)) )
      return( "featureData does not contain a 'chr' column.")
   if( ! is( fData(object)$chr, "factor" ) )
      return( "The 'chr' column in fData is not a factor." )

   if( ! "start"  %in% names(fData(object)) )
      return( "featureData does not contain a 'start' column.")
   if( ! is( fData(object)$start, "integer" ) )
      return( "The 'start' column in fData is not integer." )

   if( ! "end"  %in% names(fData(object)) )
      return( "featureData does not contain a 'end' column.")
   if( ! is( fData(object)$end, "integer" ) )
      return( "The 'end' column in fData is not integer." )

   if( ! "strand"  %in% names(fData(object)) )
      return( "featureData does not contain a 'strand' column.")
   if( ! is( fData(object)$strand, "factor" ) )
      return( "The 'strand' column in fData is not a factor." )
#   if( any( sort( levels(fData(object)$strand)) != c( "-", "+" ) ) )
#      return( "The 'strand' column in fData does not have the levels '+' and '-'." )   ### strange reason, make pasilla check crash only in the check, not sourcing exactly the same code
   if( !is(fData(object)$dispersion, "numeric")){
      return( "The 'dispersion' is not numeric")}
   if( !is(fData(object)$dispersion_CR_est, "numeric")){
      return( "The 'dispersion_CR_est' column is not numeric")}
  if( !is(fData(object)$pval, "numeric")){
      return( "The 'pval' values are not numeric")}
   if( !is.integer( assayData(object)[["counts"]] ) )
      return( "The count data is not in integer mode." )

   if( any( assayData(object)[["counts"]] < 0 ) )
      return( "The count data contains negative values." )

   if( length( object@dispFitCoefs ) != 2 )
      return( "dispFitCoefs is not a vector of length 2." )

   TRUE
} )

counts <- function( ecs ) {
   stopifnot( is( ecs, "ExonCountSet" ) )
   assayData( ecs )[[ "counts" ]]
}

`counts<-` <- function( ecs, value ) {
   stopifnot( is( ecs, "ExonCountSet" ) )
   assayData( ecs )[[ "counts" ]] <- value
   validObject( ecs )
   ecs
}   

sizeFactors <- function( ecs ) {
   stopifnot( is( ecs, "ExonCountSet" ) )
   sf <- pData(ecs)$sizeFactor
   names( sf ) <- colnames( counts(ecs) )
   sf
}

`sizeFactors<-` <- function( ecs, value ) {
   stopifnot( is( ecs, "ExonCountSet" ) )
   pData(ecs)$sizeFactor <- value
   validObject( ecs )
   ecs
}
      
design <- function( ecs, drop=TRUE, asAnnotatedDataFrame=FALSE ) {
   stopifnot( is( ecs, "ExonCountSet" ) )
   if( asAnnotatedDataFrame )
      return( phenoData(ecs)[, ecs@designColumns ] )
   ans <- pData(ecs)[, ecs@designColumns, drop=FALSE ]
   if( ncol(ans) == 1 && drop ) {
      ans <- ans[,1]
      names(ans) <- colnames( counts(ecs) ) }
   else
      rownames( ans ) <- colnames( counts(ecs) )
   ans
}

`design<-` <- function( ecs, value ) {
   stopifnot( is( ecs, "ExonCountSet" ) )
   
   # Is it multivariate or just a vector?
   if( ncol(cbind(value)) > 1 )
      value <- as( value, "AnnotatedDataFrame" )
   else {
      value <- new( "AnnotatedDataFrame",
         data = data.frame( condition = value ) )
      varMetadata( value )[ "condition", "labelDescription" ] <-
         "experimental condition, treatment or phenotype" }

   rownames( pData(value) ) <- rownames( pData(ecs) )
   dimLabels( value ) <- dimLabels( phenoData(ecs) )
   phenoData(ecs) <- combine( 
      phenoData(ecs)[ , !( colnames(pData(ecs)) %in% ecs@designColumns ), drop=FALSE ], 
      value )
   ecs@designColumns <- colnames( pData(value) )   
   validObject(ecs)
   ecs
}
      
geneIDs <- function( ecs ) {
   stopifnot( is( ecs, "ExonCountSet" ) )
   g <- fData(ecs)$geneID
   names(g) <- rownames( counts(ecs) )
   g
}

`geneIDs<-` <- function( ecs, value ) {
   stopifnot( is( ecs, "ExonCountSet" ) )
   fData(ecs)$geneIDs <- value
   validObject(ecs)
   ecs
}
      
exonIDs <- function( ecs ) {
   stopifnot( is( ecs, "ExonCountSet" ) )
   g <- fData(ecs)$exonID
   names(g) <- rownames( counts(ecs) )
   g
}

`exonIDs<-` <- function( ecs, value ) {
   stopifnot( is( ecs, "ExonCountSet" ) )
   fData(ecs)$exonIDs <- value
   validObject(ecs)
   ecs
}
            
conditions <- design
`conditions<-` <- `design<-`            


subsetByGenes <- function( ecs, genes ) {
   stopifnot( is( ecs, "ExonCountSet" ) )
   stopifnot( all( genes %in% unique(geneIDs(ecs)) ) )
   ecs2 <- ecs[ geneIDs(ecs) %in% genes, ]
   ecs2
}
   
geneCountTable <- function( ecs ) {
   stopifnot( is( ecs, "ExonCountSet" ) )
   do.call( rbind, 
      tapply( 1:nrow(ecs), geneIDs(ecs), function(rows) 
         colSums( counts(ecs)[rows,,drop=FALSE] ) ) )
}

DEUresultTable <- function(ecs, foldChange=FALSE){

	result <- data.frame(geneID=geneIDs(ecs), 
			     exonID=exonIDs(ecs), 
			     dispersion_CR_est=featureData(ecs)$dispersion_CR_est, 
			     dispersion=featureData(ecs)$dispersion, 
			     pvalue=fData(ecs)$pvalue, 
			     padjust=fData(ecs)$padjust)

	if(length(levels(design(ecs))) == 2 && foldChange){
		result$log2change <- rep(NA, nrow(result))
		for(geneID in unique(geneIDs(ecs))){
			rt<-which(featureData(ecs)$geneID==geneID)
			if(length(rt) < 2){next}
			count <- t(t(counts(ecs)[rt,])/sizeFactors(ecs))
			numcond<-length(unique(design(ecs)))
			numexons<-nrow(count)
			mf<-modelFrameForGene(ecs, geneID)
			mf$offset <- log(mf$sizeFactor)
			fam <- negative.binomial( 1 / commonDispersion(ecs) )   # <<-- TODO Replace commonDispersion
			fam$family = "Negative Binomial(varying)"
			fit1 <- try(glm( count~condition*exon, mf, family = fam, offset = mf$offset ))
			if(inherits(fit1, "try-error")){
				next
			}else{
				intercept<-fit1$coefficients["(Intercept)"]
				treatments<-fit1$coefficients[2:numcond]
				treatments<-treatments[order(names(treatments))]
				untr<-fit1$coefficients[(numcond+1):(numexons+numcond-1)]
				treatedcoef<-fit1$coefficients[(numcond+numexons):length(fit1$coefficients)]
				treatedcoef<-treatedcoef[order(sapply(strsplit(names(treatedcoef), ":"), "[[", 1))]
				coeff<-data.frame(c(intercept, untr+intercept))
				ini<-seq(1, 1+((numexons-1)*(numcond-1)), numexons-1)
				end<-seq(numexons-1, ((numexons-1)*(numcond-1)), numexons-1)
				for(i in 1:length(end)){
				coeff<-cbind(coeff, c(treatments[i]+intercept, treatedcoef[ini[i]:end[i]]+intercept+treatments[i]+untr))
				}
				result$log2change[result$geneID %in% geneID] <- log2(exp(coeff[,2])/exp(coeff[,1]))
			}
		}
	}
	result
}
