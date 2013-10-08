rmDepCols <- function(m) {
   q <- qr( m )
   if( q$rank < ncol(m) )
      m[ , -q$pivot[ (q$rank+1) : ncol(m) ] ]
   else
      m
}

constructModelFrame <- function( ecs ){
  stopifnot( inherits( ecs, "ExonCountSet" ) )
  modelFrame <- cbind(
     sample = sampleNames(ecs),
     design(ecs, drop=FALSE),
     sizeFactor = sizeFactors(ecs) )
  modelFrame <- rbind(
     cbind( modelFrame, exon="this" ),
     cbind( modelFrame, exon="others" ) )
  rownames(modelFrame) <- NULL
  return( modelFrame )
}

getCountVector <- function( ecs, geneID, exonID ) {
   stopifnot( inherits( ecs, "ExonCountSet" ) )
   ctfg <- countTableForGene( ecs, geneID )
   w <- which( rownames(ctfg) == exonID )
   stopifnot( length(w) == 1 )
   c( ctfg[ w, ], colSums( ctfg[ -w, , drop=FALSE ] ) )
}

estimateExonDispersion <- function( ecs, geneID, exonID, modelFrame, mm ){
   stopifnot( inherits( ecs, "ExonCountSet" ) )
   if( all( is.na( sizeFactors( ecs ) )) ){
     stop("Please calculate size factors before estimating dispersions\n")
   }
   
   count <- getCountVector( ecs, geneID, exonID )
   disp <- .1
   for( i in 1:10 ) {
     fit <- glmnb.fit( mm, count, disp, log( modelFrame$sizeFactor ) )
     olddisp <- disp
     disp <- exp( optimize( function(logalpha)
        profileLogLikelihood( exp(logalpha), mm, count, fitted.values(fit) ), 
        log( c( 1e-11, 1e5 ) ), maximum=TRUE, tol=.01 )$maximum )
     if( abs( log(disp) - log(olddisp) ) < .03 )
        break
     }
  disp
}

setMethod("estimateDispersions", signature(object="ExonCountSet"),
function( object, formula= ~ sample + exon + condition : exon, minCount=10, nCores=1 ){
  stopifnot( inherits( object, "ExonCountSet" ) )
   if( all( is.na( sizeFactors( object )) ) ){
     stop("Please calculate size factors before estimating dispersions\n")
   }
   testable <- rowSums(counts(object)) >= minCount
   for( r in split( seq_len(nrow(object)), geneIDs(object) ) ) {
     if( sum( testable[r] ) <= 1 )
       testable[r] <- FALSE
   }
   fData(object)$testable <- testable

   if( nCores > 1 ){
      if (!is.loaded("mc_fork", PACKAGE = "parallel")) {
          stop("Please load first parallel package or set parameter nCores to 1.")
      }else{
          myApply <- function(X, FUN){ parallel::mclapply( X, FUN, mc.cores=nCores ) }
      }
   }else{
     myApply <- lapply
   }
   object@formulas[["formulaDispersion"]] <- deparse(formula)

   rows <- seq_len(nrow(object))
   modelFrame <- constructModelFrame( object )
   mm <- rmDepCols( model.matrix( formula, modelFrame ) )
   disps <- myApply( rows, function(i) {
      if( i %% 100 == 0 )
         cat(".")
      if( fData(object)$testable[i] ) {
         a <- try( estimateExonDispersion( object, geneIDs(object)[i], exonIDs(object)[i], modelFrame, mm ), silent=TRUE)
         if( inherits( a, "try-error" ) ) {
            warning( sprintf("Unable to estimate dispersions for %s:%s", as.character( geneIDs(object)[i] ), exonIDs(object)[i]) )
            NA }
         else{
            a }
      }else{
          NA 
      }
    })
    names(disps) <- featureNames(object)
    fData(object)[names(disps), "dispBeforeSharing"] <- unlist(disps)
    fData(object)$testable[which( is.na( fData(object)$dispBeforeSharing ) )] <- FALSE
    cat("\nDone\n")
    object
} )

testExonForDEU <- function(ecs, geneID, exonID, modelFrame, mm0, mm1, disp){
  stopifnot( inherits( ecs, "ExonCountSet" ) )
  stopifnot( !is.na(disp) )
  stopifnot( any(geneIDs(ecs) %in% geneID & exonIDs(ecs) %in% exonID ) )
  if( all( is.na( sizeFactors( ecs )) ) ){
    stop("Please calculate size factors first\n")
  }
  count <- getCountVector( ecs, geneID, exonID )
  fit0  <- glmnb.fit( mm0,  count, disp, log( modelFrame$sizeFactor ) )
  fit1 <- glmnb.fit( mm1, count, disp, log( modelFrame$sizeFactor ) )
  pchisq( deviance( fit0 ) - deviance( fit1 ), ncol( mm1 ) - ncol( mm0 ), lower.tail=FALSE )
}

testForDEU <- function( ecs, formula0 = ~ sample + exon, 
    formula1 = ~ sample + exon + condition:exon, dispColumn="dispersion", nCores=1 ){
  stopifnot( inherits( ecs, "ExonCountSet" ) )
   if( all( is.na( sizeFactors( ecs )))) {
     stop("Please calculate size factors before estimating dispersions\n")
   } 
   if( all( is.na( fData(ecs)[,dispColumn] ) ) ){
     stop("Please estimate dispersions before calling this function\n")
   }
   if( nCores > 1 ){
      if (!is.loaded("mc_fork", PACKAGE = "parallel")) {
          stop("Please load first parallel package or set parameter nCores to 1...")
      }else{
          myApply <- function(X, FUN){ parallel::mclapply( X, FUN, mc.cores=nCores ) }
      }
   }else{
     myApply <- lapply
   }
   ecs@formulas[["formula0"]] <- deparse(formula0)
   ecs@formulas[["formula1"]] <- deparse(formula1)
   rows <- seq_len(nrow(ecs))
   modelFrame <- constructModelFrame( ecs )
   mm0 <- rmDepCols( model.matrix( formula0, modelFrame ) )
   mm1 <- rmDepCols( model.matrix( formula1, modelFrame ) )

   pvals <- myApply( rows,
     function(i) {
       if( i %% 1000 == 0 )
          cat(".")
       if( fData(ecs)$testable[i] ) {
          a <- try( testExonForDEU( ecs, geneIDs(ecs)[i], exonIDs(ecs)[i], modelFrame, mm0, mm1, fData(ecs)[i, dispColumn] ) )
          if( any(inherits( a, "try-error" ) )) {
             warning( sprintf("Unable to calculate p-values for %s:%s\n", as.character( geneIDs(ecs)[i] ), exonIDs(ecs)[i]) )
             NA }
          else
             a }
       else
          NA })

    names( pvals ) <- featureNames( ecs )
    fData(ecs)[names(pvals), "pvalue"] <- unlist(pvals)
    fData(ecs)$padjust <- p.adjust( fData(ecs)$pvalue, method="BH" )
    ecs

}




doCompleteDEUAnalysis <- function( ecs, formula0 = ~ sample + exon, formula1 = ~ sample + exon + condition : exon, minCount=10,
     nCores=1, path=NULL, FDR=0.1, fitExpToVar="condition", color=NULL, color.samples=NULL )
{
   stopifnot(is(ecs, "ExonCountSet"))
   ecs <- estimateSizeFactors( ecs )
   ecs <- estimateDispersions( ecs, formula1, nCores=nCores, minCount=minCount )
   ecs <- fitDispersionFunction( ecs )
   ecs <- testForDEU( ecs, formula0=formula0, formula1=formula1, nCores=nCores )
   ecs <- estimatelog2FoldChanges(ecs, fitExpToVar=fitExpToVar, nCores=nCores )
   if( !is.null(path) ) {
      DEXSeqHTML( ecs, path=path, FDR=FDR, fitExpToVar=fitExpToVar, color=color, color.samples=color.samples )
   }
   ecs
}

