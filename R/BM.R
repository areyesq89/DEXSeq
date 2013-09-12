estimateExonDispersionsForModelFrame_BM <- function( modelFrame, formula=NULL, mm=NULL, muhat=NULL, initialGuess = .01 )
{
   exonNames <- as.character(levels(modelFrame$exon))
   if( length( exonNames ) <= 1 ) {
      ans <- NA_real_
      names(ans) <- exonNames
      return( ans )
   }
   if(is.null(formula)){
	formula <- count ~ sample + condition*exon
   }

   if(is.null(mm)){
      mm <- model.matrix( formula, modelFrame )

   }


   countsums <- tapply( modelFrame$count, modelFrame$exon, sum )
   y <- modelFrame$count

   mm <- rmDepCols( mm )

   if(is.null(muhat)){
      if( nrow(mm) <= ncol(mm) )
         stop( "Underdetermined model; cannot estimate dispersions. Maybe replicates have not been properly specified." )
      muhat <- fitted.values(glmnb.fit(mm, y, initialGuess, log(modelFrame$sizeFactor)))
   }

   disp <- rep( initialGuess, length(exonNames) )
   names(disp) <- exonNames
      for( exon in exonNames[ countsums > 0 ] )
         disp[exon] <- exp(
            optimize(
               function(logalpha) {
                  disp[exon] <- exp( logalpha )
                  profileLogLikelihood( disp[as.character(modelFrame$exon)], mm, y, muhat ) },
               log( c( 1e-11, 1e5 ) ),
   	       tol = 0.1,
               maximum=TRUE
            )$maximum )
   disp[ countsums == 0 ] <- NA
   disp[ disp < 1e-10 ] <- 0
   disp
}


estimateDispersions_BM <- function( object, formula=count ~ sample + condition*exon, 
	   initialGuess=.01, nCores=1, minCount=10, maxExon=70, quiet=FALSE, file="")
   {
      cds <- object
      stopifnot(is(cds, "ExonCountSet"))
      if( all( is.na(sizeFactors(cds)) ) ){
         stop( "Estimate size factors first." )
      }
      cds@formulas[["formulaDispersion_BM"]] <- deparse(formula)


      if(!interactive() & !quiet & file == ""){  ## if the session is not interactive, quiet is FALSE and there is no file, then quiet does not make any sense
         quiet=TRUE
      }
      ########## DEFINE TESTABLE GENES #############
      # Which exons are actually testable?
       # take away those exons with counts lower than minCounts
      testable <- rowSums(counts(cds)) >= minCount
      if(!all(testable) & nCores<=1){
         warning(sprintf("Exons with less than %d counts will not be tested. For more details please see the manual page of 'estimateDispersions', parameter 'minCount'", minCount))
      }
      # If a gene contains less than two high count exons, all its exons non-testable
      for( r in split( seq_len(nrow(cds)), geneIDs(cds) ) ) {
         if( sum( testable[r] ) <= 1 | sum( testable[r] ) > maxExon )
            testable[r] <- FALSE
      }
      fData(cds)$testable <- testable

      # take away genes with just one exon
      generle <- rle( sort(as.character( geneIDs(cds) ) ))
#      fData(cds)$testable[which(geneIDs(cds) %in% generle$values[which( generle$lengths ==1 )])] <- FALSE
      # take away those exons bigger than maxExon (default 70)
#      fData(cds)$testable[which(geneIDs(cds) %in% generle$values[which( generle$lengths > maxExon )])] <- FALSE
      ###

      if(max(generle$lengths) > maxExon & nCores<=1) {
         warning(sprintf("Genes with more than %d testable exons will be omitted from the analysis. For more details please see the manual page of 'estimateDispersions', parameter 'maxExon'.", maxExon))
      }

      ##### DOES THE SAME AS A TAPPLY, but without considering the order of the levels ##
#      testable <- fData(cds)$testable
      testablegenes <- as.character(unique(fData(cds)[which(fData(cds)$testable),]$geneID))

      if(!quiet & nCores==1 ) {
         cat( "Dispersion estimation. (Progress report: one dot per 100 genes)\n", file=file, append=TRUE)
      }

      i <- 0
      if(nCores > 1){
         if(!quiet){
         cat(sprintf("Estimating Cox-Reid exon dispersion estimates using %d cores. (Progress report: one dot per 100 genes)\n", nCores), file=file, append=TRUE)}
         toapply <- function(x){estimateDispersions(x, formula=formula, initialGuess=initialGuess, nCores=-1, minCount=minCount, maxExon=maxExon, file=file, quiet=quiet)}
         cds <- divideWork(cds, funtoapply=toapply, fattr="dispBeforeSharing", mc.cores=nCores, testablegenes)
      }else{
         modelFrames <- lapply( testablegenes, function(gs) {
            mf <- modelFrameForGene(cds, gs, onlyTestable=TRUE)
            mf$offset <- log(mf$sizeFactor)
            mf })

         names(modelFrames) <- testablegenes

         modelmatrices <- lapply( modelFrames, function(mf){
            model.matrix( formula, data = mf )} )

         names(modelmatrices) <- testablegenes

         muhats <- lapply( testablegenes, function( gn ) {
            y <- modelFrames[[gn]]$count
            y1 <- pmax(y, 1/6)
            mf <- modelFrames[[gn]]
            mm <- modelmatrices[[gn]]
            mm <- rmDepCols(mm)
            if( nrow(mm) <= ncol(mm) )
               stop( "Underdetermined model; cannot estimate dispersions. Maybe replicates have not been properly specified." )
            muhat <- try(fitted.values(glmnb.fit(mm, y, initialGuess, mf$offset)), silent=TRUE)
            muhat
         })

         names(muhats) <- testablegenes

         badones <- which( sapply( muhats, inherits, "try-error") )
         if( length(badones) > 0 ) {
            testablegenes <- testablegenes[ ! testablegenes %in% names(badones) ]
    	      warning( paste( "Failed to set up model frames for genes ",
   	         paste( names(badones), collapse=", " ) ) )
         }

         ###### WORKS FASTER THIS WAY THAN IN EVERY CYCLE ACCESSING fData
         fData(cds)$dispBeforeSharing <- NA_real_
         dispBeforeSharing <- fData(cds)$dispBeforeSharing
         testable <- fData(cds)$testable
         exonids <- fData(cds)$exonID

         for( genename in testablegenes ){
            mf <- modelFrames[[genename]]
            disps <- try(
               estimateExonDispersionsForModelFrame_BM( mf,
                  mm=modelmatrices[[genename]], muhat=muhats[[genename]] ),
               silent=TRUE )
            rows <- as.character(geneIDs(cds)) %in% genename & testable
            if( inherits( disps, "try-error" ) ) {
               disps <- rep( NA_real_, sum( rows ) )
               warning( sprintf( "Failed to fit dispersion for gene %s", genename ) )
            }
            stopifnot(all(names(disps)==exonids[rows]))
            dispBeforeSharing[rows] <- disps

            i <- i + 1
            if(!quiet & i %% 100 == 0 ){
               cat( ".", file=file, append=TRUE)}
         }
         fData(cds)$dispBeforeSharing <- dispBeforeSharing
      }
      if(nCores==1 & !quiet){
         cat("\n", file=file, append=TRUE)
      }
      cds
   }


testGeneForDEU_BM <- function (ecs, gene, formula0=NULL, formula1=NULL ){
   stopifnot(is(ecs, "ExonCountSet"))
   if( all( is.na(featureData(ecs)$dispersion ) ) ) {
      stop("No dispersion values found, call function fitDispersionFunction first.")
   }

   mf <- modelFrameForGene(ecs, gene, onlyTestable=TRUE)

   ans <- data.frame( row.names = levels( mf$exon ) )
   ans$deviance = NA_real_
   ans$df = NA_integer_
   ans$pvalue = NA_real_

   if(is.null(formula0)){
      formula0 <- count ~ sample + exon + condition}
   if(is.null(formula1)){
      formula1 <- count ~ sample + exon + condition * I(exon == exonID)}


   # This makes sure that the formula see the 'exonID' variable used
   # below even if it the formula was supplied as parameter
   environment(formula0) <- environment()
   environment(formula1) <- environment()

   mm <- model.matrix(formula0, mf)
   y1 <- pmax(mf$count, 1/6)
   fit <- lm.fit(mm, log(y1) - log(mf$sizeFactor))
   mm <- mm[,!is.na(fit$coefficients)]
   start <- fit$coefficients[!is.na(fit$coefficients)]

   fit0 <- try(
      glmnb.fit(mm, mf$count, mf$dispersion, log(mf$sizeFactor), coef.start=start),
   silent=TRUE)
   if( inherits( fit0, "try-error" ) ) {
      warning( sprintf( "Error in fit0 for gene %s: %s", gene, fit0 ) )
      return(ans) }

   for( exonID in unique(as.character(mf$exon)) ) {
      mm <- model.matrix(formula1, mf)
      y1 <- pmax(mf$count, 1/6)
      fit <- lm.fit(mm, log(y1) - log(mf$sizeFactor))
      mm <- mm[,!is.na(fit$coefficients)]
      start <- fit$coefficients[!is.na(fit$coefficients)]

      fit1 <- try(
         glmnb.fit(mm, mf$count, mf$dispersion, log(mf$sizeFactor), coef.start=start),
      silent=TRUE)
      if( inherits( fit1, "try-error" ) ) {
         warning( sprintf( "Error in fit1 for gene %s, exon %s: %s", gene, exonID, fit1 ) )
         next }
      ans[ exonID, "deviance" ] <- deviance( fit0 ) - deviance( fit1 )
      ans[ exonID, "df" ] <- length(fit1$coefficients) - length(fit0$coefficients)
      ans[ exonID, "pvalue"   ] <- 1 - pchisq( ans[exonID,"deviance"], ans[exonID,"df"] )
  }
  ans
}


testForDEU_BM <- function( ecs, formula0=NULL, formula1=NULL, nCores=1, quiet=FALSE, file="")
{
   stopifnot(is(ecs, "ExonCountSet"))
   if( all( is.na(featureData(ecs)$dispersion ) ) ) {
      stop("No dispersion values found, call function fitDispersionFunction first.")
   }

   if(!interactive() & !quiet & file == ""){  ## if the session is not interactive, quiet is FALSE and there is no file, then quiet does not make any sense
      quiet=TRUE
   }

   if(!quiet & nCores==1){
      cat( "Testing for differential exon usage. (Progress report: one dot per 100 genes)\n", file=file, append=TRUE)}

   testablegenes <- as.character(unique(fData(ecs)[which(fData(ecs)$testable),]$geneID))

   if(nCores > 1){
      if(!quiet){
      cat(sprintf("Testing for differential exon usage using %d cores. (Progress report: one dot per 100 genes)\n", nCores), file=file, append=TRUE)}
      toapply <- function(x){testForDEU_BM(x, formula0=formula0, formula1=formula1, nCores=-1, file=file, quiet=quiet)}
      ecs <- divideWork(ecs, funtoapply=toapply, fattr="pvalue", mc.cores=nCores, testablegenes)
   }else{
      i <- 0
      testable <- fData(ecs)$testable
      pvalue <- fData(ecs)$pvalue
      exonids <- fData(ecs)$exonID
      for( genename in testablegenes ) {
         i <- i + 1
         if(!quiet & i %% 100 == 0 ){
            cat( ".", file=file, append=TRUE)}
         rows <- as.character(geneIDs(ecs)) %in% genename & testable
         res <- testGeneForDEU_BM( ecs, genename, formula0, formula1 )
         stopifnot(all(rownames(res) == exonids[rows]))    ### makes sure to do the assignment correctly
         pvalue[rows] <- res$pvalue
      }
      fData(ecs)$pvalue <- pvalue
   }
   if(nCores==1 & !quiet){
      cat("\n", file=file, append=TRUE)
   }
   fData(ecs)$padjust <- p.adjust(fData(ecs)$pvalue, method="BH")
   ######### STORE FORMULAS IN THE OBJECTs
   if(is.null(formula0)){
      formula0 <- count~sample+exon+condition
   }
   if(is.null(formula1)){
      formula1 <- count~sample+exon+condition*I(exon==exonID)
   }
   ecs@formulas[["formula0_BM"]] <- deparse(formula0)
   ecs@formulas[["formula1_BM"]] <- deparse(formula1)
   ecs
}

makeCompleteDEUAnalysis_BM <- function(ecs, formulaDispersion=count ~ sample + condition*exon, minCount=10, maxExon=50, formula0=NULL, formula1=NULL, FDR=0.1, fitExpToVar="condition", nCores=1, path=NULL, color=NULL, color.samples=NULL, quiet=FALSE, file="")
{
   stopifnot(is(ecs, "ExonCountSet"))
   ecs <- estimateSizeFactors( ecs )
   ecs <- estimateDispersions_BM( ecs, formulaDispersion, nCores=nCores, quiet=quiet, file=file, minCount=minCount, maxExon=maxExon)
   ecs <- fitDispersionFunction( ecs )
   ecs <- testForDEU_BM( ecs, formula1=formula1, formula0=formula0, nCores=nCores, quiet=quiet, file=file)
   ecs <- estimatelog2FoldChanges(ecs, fitExpToVar=fitExpToVar, nCores=nCores, quiet=quiet, file=file)
   if(!is.null(path)){
      DEXSeqHTML(ecs, path=path, FDR=FDR, fitExpToVar=fitExpToVar, color=color, color.samples=color.samples)
   }
   ecs
}
