estimateSizeFactors <- function( ecs ){
   stopifnot( is( ecs, "ExonCountSet") )
   geomeans <- exp( rowMeans( log( counts(ecs) ) ) )
   sizeFactors(ecs) <- 
      apply( counts(ecs), 2, function(cnts) 
         median( ( cnts / geomeans )[ geomeans>0 ] ) )
   ecs
}

countTableForGene <- function( ecs, geneID, normalized=FALSE, withDispersion=FALSE ) {
   stopifnot( is( ecs, "ExonCountSet") )
   rows <- geneIDs(ecs) == geneID
   if( sum(rows) == 0 )
      stop( "geneID not found" )
	ans <- counts( ecs )[ rows, , drop=FALSE ]
	rownames( ans ) <- exonIDs(ecs)[ rows ]
	attr( ans, "geneID" ) <- geneID
	if( normalized )
	   if(any(is.na(sizeFactors(ecs)))){
   	   stop("Please first call function estimateSizeFactors\n")
	   }else{
	   ans <- t( t(ans) / sizeFactors(ecs) )
	   }
	if( withDispersion )
      attr( ans, "dispersion" ) <- fData(ecs)$dispersion[rows]
   ans
}

modelFrameForGene <- function( ecs, geneID ) {
   stopifnot( is( ecs, "ExonCountSet") )
   rows <- geneIDs(ecs) == geneID
   numExons <- sum(rows)
   exonCol <- rep( factor( exonIDs(ecs)[rows] ), ncol( counts(ecs) ) )
   modelFrame <- data.frame(
      sample = rep( factor( colnames( counts(ecs) ) ), each = numExons ),
      exon = exonCol,
      sizeFactor = rep( sizeFactors(ecs), each = numExons ) )
   for( cn in colnames( design(ecs,drop=FALSE) ) )
      modelFrame[[cn]] <- rep( design(ecs,drop=FALSE)[[cn]], each=numExons )
   modelFrame$dispersion <- fData(ecs)$dispersion[ rows ][ 
      match( modelFrame$exon, exonIDs(ecs)[rows] ) ]
   modelFrame$count <- as.vector( counts(ecs)[rows,] )
   attr( modelFrame, "geneID" ) <- geneID
   modelFrame
}

assertOneWay <- function( ecs ) {
   stopifnot( is( ecs, "ExonCountSet") )
   if( ncol( design(ecs,drop=FALSE) ) != 1 )
      stop( "This function only works for simple one-way designs, but your have specified a design with multiple predictors" )
   if( colnames( design(ecs,drop=FALSE) ) != "condition" )
      stop( "This function requires the design column to have its default name, 'condition'." )
}


profileLogLikelihood <- function( disp, mm, y, muhat )
{
   # calculate the log likelihood:
   if(length(disp) != length(y)){
      disp <- rep(disp, length(y))
   }

   ll <- sum( sapply( 1:length(y), function(i) 
      dnbinom( y[i], mu=muhat[i], size=1/disp[i], log=TRUE ) ) )

   # transform the residuals, i.e., y - muhat, to the linear
   # predictor scale by multiplying them with the derivative
   # of the link function, i.e., by 1/muhat, and add this to the
   # linear predictors, log(muhat), to get the predictors that
   # are used in the IWLS regression
   z <- log(muhat) + ( y - muhat ) / muhat

   # the variance function of the NB is as follows
   v0 <- muhat + disp * muhat^2 

   # transform the variance vector to linear predictor scale by
   # multiplying with the squared derivative of the link to
   # get the (reciprocal) weights for the IWLS
   w <- 1 / ( ( 1 / muhat )^2 * v0 )

   # All we need from the IRLS run is the QR decomposition of
   # its matrix
   qrres <- qr( mm*sqrt(w) )

   # from it, we extract we leverages and calculate the Cox-Reid
   # term:
   cr <- sum( log( abs( diag( qrres$qr )[ 1:qrres$rank ] ) ) )  

   # return the profile log likelihood:
   ll - cr 
}

estimateExonDispersionsForModelFrame <- function( modelFrame, formula=NULL, mm=NULL, muhat=NULL, initialGuess = .01 )
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

   if(is.null(muhat)){
      muhat <- fitted.values( glm.fit( mm, y, family=negative.binomial(1/initialGuess), offset=log(modelFrame$sizeFactor) ) )
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
               maximum=TRUE 
            )$maximum )
   disp[ countsums == 0 ] <- NA
   disp[ disp < 1e-10 ] <- 0
   disp
}

fitDispersions <- function( ecs )
{
   means <- colMeans( t(counts(ecs))/sizeFactors(ecs) )
   disps <- fData(ecs)$dispersion_CR_est
   coefs <- c( .1, 1 )
   iter <- 0   
   while(TRUE) {
      residuals <- disps / ( coefs[1] + coefs[2] / means )
      good <- (residuals > 1e-4) & (residuals < 15)
      fit <- glm( disps[good] ~ I(1/means[good]), 
         family=Gamma(link="identity"), start=coefs )
      oldcoefs <- coefs   
      coefs <- coefficients(fit)
      if( sum( log( coefs / oldcoefs )^2 ) < .005 )
         break
      iter <- iter + 1
      if( iter > 10 ) {
         warning( "Dispersion fit did not converge." )
         break }
    }
  return(coefs)
}


estimateDispersions <- function( ecs, formula=NULL, file=NULL, initialGuess=0.01, quiet=FALSE)
{
   stopifnot(is(ecs, "ExonCountSet"))
   if( all( is.na(sizeFactors(ecs)) ) ){
      stop( "Estimate size factors first." )	
   }
   
   # Which exons are actually testable? 
   # An all-zero gene is not testable
   fData(ecs)$testable <- rowSums(counts(ecs)) > 0
   # If a gene contains less then two non-zero exons, all its exons non-testable
   fData(ecs)$testable <- unlist( tapply( fData(ecs)$testable, geneIDs(ecs), function(x) 
      if( sum(x) > 1 ) x else rep( FALSE, length(x) ) ) ) 
   testableGenes <- names( which( tapply( fData(ecs)$testable, geneIDs(ecs), any ) ) )

   if( !quiet ) {
      cat( "Dispersion estimation. (Progress report: one dot per 100 genes)\n" )
      cat( "  Setting up model frames ")
   }
   i <- 0
   modelFrames <- sapply( testableGenes, function(gs) {
      mf <- modelFrameForGene(ecs, gs)
      mf$offset <- log(mf$sizeFactor)
      if( !quiet ) {
         i <<- i + 1
         if( i %% 100 == 0 )
            cat( "." )
      }
      mf }, simplify=FALSE )
   
   
   if( is.null(formula) ){
      formula <- count ~ sample + condition*exon
   }
   
   
   if( !quiet )
      cat( "\n  Setting up model matrices." )
   modelmatrices <- sapply( modelFrames, function(mf)
      model.matrix( formula, data = mf ), simplify=FALSE )

   if( !quiet )
      cat( "\n  Calculating initial fits." )

   # warnOpt <- getOption( "warn" )
   # options( warn=-1 )  # TODO: Write a proper handler to silence 
      
   muhats <- sapply( testableGenes, function( gn ) { 
      try( 
        fitted.values(glm.fit( 
            modelmatrices[[gn]], modelFrames[[gn]]$count, 
            family = negative.binomial(1/initialGuess), 
            offset = modelFrames[[gn]]$offset )),
         silent=TRUE ) },
      simplify=FALSE )
   # options( warn=warnOpt )

   badones <- which( sapply( muhats, inherits, "try-error") )
   if( length(badones) > 0 ) {
      testableGenes <- testableGenes[ ! testableGenes %in% names(badones) ]
	   warning( paste( "Failed to set up model frames for genes ", 
	      paste( names(badones), collapse=", " ) ) )
   }
   
   fData(ecs)$dispersion_CR_est <- NA_real_
   
   if( !quiet ) {
      cat( "\n  Performing Cox-Reid dispersion estimation " )
      i <- 0 }
   for( genename in testableGenes ){

      mf <- modelFrames[[genename]]
      # warnOpt <- getOption( "warn" )
      # options( warn=-1 )  # TODO: Write a proper handler to silence warnings
      disps <- try( 
         estimateExonDispersionsForModelFrame( mf, 
            mm=modelmatrices[[genename]], muhat=muhats[[genename]] ),
         silent=TRUE )
      # options( warn=warnOpt )
      if( inherits( disps, "try-error" ) ) {
         disps <- rep( NA_real_, length( muhats[[genename]] ) )
         warning( sprintf( "Failed to fit dispersion for gene %s", genename ) )
      }   

      if( !is.null(file) ) {
         countsums <- tapply( mf$count, mf$exon, sum )
         stopifnot( identical( names(disps), names(countsums) ) )   # What is this for?
         write.table( 
            data.frame( gene=genename, exon=names(disps), disp=disps, countsum=countsums ), 
            quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t", file=file )
         flush( file ) }

      rows <- as.character(geneIDs(ecs)) %in% genename
#     stopifnot( identical( names(disps), unname(exonIDs(ecs)[rows])) ) for the momment
      fData(ecs)$dispersion_CR_est[rows] <- disps 
      
      if( !quiet ) {
         i <- i + 1
         if( i %% 100 == 0 )
            cat( "." ) }      
   }

   if( !quiet )
      cat( "\n  Fitting mean-dispersion relation" )
   ecs@dispFitCoefs <- fitDispersions( ecs )  
   
   fData(ecs)$dispersion <- pmin(
      pmax( 
         fData(ecs)$dispersion_CR_est, 
         ecs@dispFitCoefs[1] + ecs@dispFitCoefs[2] / colMeans( t(counts(ecs))/sizeFactors(ecs) ),
         na.rm = TRUE ),
      1e8 )   # 1e8 as an arbitrary way-too-large value to capture infinities
 
   if( !quiet )
      cat( "\nFinished with dispersion estimation.\n" )
   ecs
}


testGeneForDEU <- function (ecs, gene, formula0=NULL, formula1=NULL, 
      glm.control = list( maxit=100, epsilon=3e-4 ) ){
   stopifnot(is(ecs, "ExonCountSet"))
   if( all( is.na(featureData(ecs)$dispersion ) ) ) {
      stop("No dispersion values found, call function estimateDispersions first.")		
   }

   mf <- modelFrameForGene(ecs, gene)
   
   ans <- data.frame( row.names = levels( mf$exon ) )
   ans$deviance = NA_real_ 
   ans$df = NA_integer_
   ans$pvalue = NA_real_

   # Which exons have no counts in any sample?
   nonzero <- tapply( mf$count, mf$exon, sum ) > 0
   
   # We need at least two non-zero exons to do a test
   if( sum(nonzero) <= 1 )
      return( ans )
      
   # Remove all the zero exons from the model frame
   mf <- mf[ nonzero[mf$exon], ]
   
   if(is.null(formula0)){
      formula0 <- count ~ sample + exon + condition
   }
   
   if(is.null(formula1)){
      formula1 <- count ~ sample + exon + condition * I(exon == exonID)
   }
   
   # This makes sure that the formula see the 'exonID' variable used
   # below even if it the formula was supplied as parameter
   environment(formula0) <- environment()
   environment(formula1) <- environment()

   for( exonID in unique(as.character(mf$exon)) ) {
      if( !nonzero[exonID] ) { 
         next }
      fam <- negative.binomial( 1 / mf$dispersion )
      fam$family = "Negative Binomial(varying)"
      fit0 <- try( 
         glm( formula0, data=mf, family = fam, offset = log(mf$sizeFactor), control=glm.control ),
         silent=TRUE )
      if( inherits( fit0, "try-error" ) ) {
         warning( sprintf( "Error in fit0 for gene %s, exon %s: %s", gene, exonID, fit0 ) )
         next }
      fit1 <- try( 
         glm( formula1, data=mf, family = fam, offset = log(mf$sizeFactor), control=glm.control ),
         silent=TRUE )
      if( inherits( fit1, "try-error" ) ) {
         warning( sprintf( "Error in fit1 for gene %s, exon %s: %s", gene, exonID, fit1 ) )
         next }
      ans[ exonID, "deviance" ] <- deviance( fit0 ) - deviance( fit1 )
      ans[ exonID, "df" ]       <- fit0$df.residual - fit1$df.residual
      ans[ exonID, "pvalue"   ] <- 1 - pchisq( ans[exonID,"deviance"], ans[exonID,"df"] )
  }
  ans
}

testForDEU <- function( ecs, formula0=NULL, formula1=NULL, padjust=TRUE, quiet=FALSE )
{
   stopifnot(is(ecs, "ExonCountSet"))
   i <- 0
   if( !quiet )
      cat( "Testing for differential exon usage " )
   for( genename in as.character( unique( geneIDs(ecs) ) ) ) {
      i <- i + 1
      if( !quiet & i %% 100 == 0 )
         cat( "." )
      rows <- as.character(geneIDs(ecs)) %in% genename
      res <- testGeneForDEU( ecs, genename, formula0, formula1 )
#      stopifnot( identical( rownames(res), unname(exonIDs(ecs)[rows])) )  ### temporal for HTSEq order of the exons
#      fData(ecs)$deviance[rows] <- res$deviance  ##not used at all
 #     fData(ecs)$df[rows] <- res$df  # not use at all
      fData(ecs)$pvalue[rows] <- res$pvalue }
   if(padjust){
      fData(ecs)$padjust <- p.adjust(fData(ecs)$pvalue, method="BH")
   }
   if( !quiet )
      cat( "\n" )
   ecs   
}

makeExampleExonCountSet <- function(){
   countfiles <- dir(system.file("files/", package="DEXSeq"), pattern="fbsubset.txt")[c(2, 3, 6, 7)]
   aggregatefile <- dir(system.file("files/", package="DEXSeq"), pattern="subset.gff")
   aggregatefile <- paste(system.file("files/", package="DEXSeq"), aggregatefile, sep="")
   countfiles <- paste(system.file("files/", package="DEXSeq"), countfiles, sep="")
   ecs <- read.HTSeqCounts(countfiles, c("treated", "treated", "untreated", "untreated"), aggregatefile=aggregatefile)
   ecs
}

