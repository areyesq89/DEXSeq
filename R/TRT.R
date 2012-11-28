rmDepCols <- function(m) {
   q <- qr( m )
   if( q$rank < ncol(m) )
      m[ , -q$pivot[ (q$rank+1) : ncol(m) ] ]
   else
      m
}

modelFrameForTRT <- function( ecs ){
  stopifnot( inherits( ecs, "ExonCountSet" ) )
  modelFrame <- cbind(
     sample = sampleNames(ecs),
     design(ecs),
     sizeFactor = sizeFactors(ecs) )
  modelFrame <- rbind(
     cbind( modelFrame, exon="this" ),
     cbind( modelFrame, exon="others" ) )
  rownames(modelFrame) <- NULL
  return( modelFrame )
}

getCountVectorTRT <- function( ecs, geneID, exonID ) {
   stopifnot( inherits( ecs, "ExonCountSet" ) )
   ctfg <- countTableForGene( ecs, geneID )
   w <- which( rownames(ctfg) == exonID )
   stopifnot( length(w) == 1 )
   c( ctfg[ w, ], colSums( ctfg[ -w, , drop=FALSE ] ) )
}

estimateExonDispersionTRT <- function( ecs, geneID, exonID, modelFrame, mm ){
   stopifnot( inherits( ecs, "ExonCountSet" ) )
   if( all( is.na( sizeFactors( ecs ) )) ){
     stop("Please calculate size factors before estimating dispersions\n")
   }
   
   count <- getCountVectorTRT( ecs, geneID, exonID )
   disp <- .1
   for( i in 1:10 ) {
     fit <- glmnb.fit( mm, count, disp, log( modelFrame$sizeFactor ) )
     olddisp <- disp
     disp <- exp( optimize( function(logalpha)
     profileLogLikelihood( exp(logalpha), mm, count, fitted.values(fit) ), log( c( 1e-11, 1e5 ) ), maximum=TRUE, tol=.01 )$maximum )
     if( abs( log(disp) - log(olddisp) ) < .03 )
        break
     }
  disp
}

estimateDispersionsTRT <- function( ecs, formula= ~ sample + condition * exon, nCores=1, minCount=10 ){
  stopifnot( inherits( ecs, "ExonCountSet" ) )
   if( all( is.na( sizeFactors( ecs )) ) ){
     stop("Please calculate size factors before estimating dispersions\n")
   }
   testable <- rowSums(counts(ecs)) >= minCount
   for( r in split( seq_len(nrow(ecs)), geneIDs(ecs) ) ) {
     if( sum( testable[r] ) <= 1 )
       testable[r] <- FALSE
   }
   fData(ecs)$testable <- testable

   if(!all(testable) & nCores<=1){
      warning(sprintf("Exons with less than %d counts will not be tested. For more details please see the manual page of 'estimateDispersions', parameter 'minCount'", minCount))
   }

   if( nCores > 1 ){
      if (!is.loaded("mc_fork", PACKAGE = "parallel")) {
          stop("Please load first parallel package or set parameter nCores to 1...")
      }else{
          myApply <- function(X, FUN){ parallel:::mclapply( X, FUN, mc.cores=nCores ) }
      }
   }else{
     myApply <- lapply
   }
   ecs@formulas[["formulaDispersionTRT"]] <- deparse(formula)

   rows <- seq_len(nrow(ecs))
   modelFrame <- modelFrameForTRT( ecs )
   mm <- rmDepCols( model.matrix( formula, modelFrame ) )
   disps <- myApply( rows, function(i) {
      if( i %% 100 == 0 )
         cat(".")
      if( fData(ecs)$testable[i] ) {
         a <- try( estimateExonDispersionTRT( ecs, geneIDs(ecs)[i], exonIDs(ecs)[i], modelFrame, mm ) )
         if( inherits( a, "try-error" ) ) {
            warning( sprintf("Unable to estimate dispersions for %s:%s", as.character( geneIDs(ecs)[i] ), exonIDs(ecs)[i]) )
            NA }
         else{
            a }
      }else{
          NA 
      }
    })
    names(disps) <- featureNames(ecs)
    fData(ecs)[names(disps), "dispBeforeSharing"] <- unlist(disps)
    fData(ecs)$testable[which( is.na( fData(ecs)$dispBeforeSharing ) )] <- FALSE
    cat("\nDone\n")
    ecs
}


testExonForDEUTRT <- function(ecs, geneID, exonID, modelFrame, mm0, mm1, disp){
  stopifnot( inherits( ecs, "ExonCountSet" ) )
  stopifnot( !is.na(disp) )
  stopifnot( any(geneIDs(ecs) %in% geneID & exonIDs(ecs) %in% exonID ) )
  if( all( is.na( sizeFactors( ecs )) ) ){
    stop("Please calculate size factors first\n")
  }
  count <- getCountVectorTRT( ecs, geneID, exonID )
  fit0  <- glmnb.fit( mm0,  count, disp, log( modelFrame$sizeFactor ) )
  fit1 <- glmnb.fit( mm1, count, disp, log( modelFrame$sizeFactor ) )
  pchisq( deviance( fit0 ) - deviance( fit1 ), ncol( mm1 ) - ncol( mm0 ), lower.tail=FALSE )
}

testForDEUTRT <- function( ecs, nCores=1, formula0= ~sample + condition + exon, formula1= ~sample + condition * exon, dispColumn="dispersion"){
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
          myApply <- function(X, FUN){ parallel:::mclapply( X, FUN, mc.cores=nCores ) }
      }
   }else{
     myApply <- lapply
   }
   ecs@formulas[["formula0TRT"]] <- deparse(formula0)
   ecs@formulas[["formula1TRT"]] <- deparse(formula1)
   rows <- seq_len(nrow(ecs))
   modelFrame <- modelFrameForTRT( ecs )
   mm0 <- rmDepCols( model.matrix( formula0, modelFrame ) )
   mm1 <- rmDepCols( model.matrix( formula1, modelFrame ) )

   pvals <- myApply( rows,
     function(i) {
       if( i %% 1000 == 0 )
          cat(".")
       if( fData(ecs)$testable[i] ) {
          a <- try( testExonForDEUTRT( ecs, geneIDs(ecs)[i], exonIDs(ecs)[i], modelFrame, mm0, mm1, fData(ecs)[i, dispColumn] ) )
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



