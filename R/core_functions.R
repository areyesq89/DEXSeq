setMethod("estimateSizeFactors", signature(cds="ExonCountSet"),
   function( cds ){
      stopifnot( is( cds, "ExonCountSet") )
      geomeans <- exp( rowMeans( log( counts(cds) ) ) )
      sizeFactors(cds) <- 
         apply( counts(cds), 2, function(cnts) 
            median( ( cnts / geomeans )[ geomeans>0 ] ) )
      cds
   }
)

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
      modelFrame[[cn]] <- factor(rep( design(ecs,drop=FALSE)[[cn]], each=numExons ))
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
   	       tol = 0.1,
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
    ecs@dispFitCoefs <- coefs
    fData(ecs)$dispersion <- pmin(
       pmax( 
          fData(ecs)$dispersion_CR_est, 
          ecs@dispFitCoefs[1] + ecs@dispFitCoefs[2] / colMeans( t(counts(ecs))/sizeFactors(ecs) ),
          na.rm = TRUE ),
          1e8 )   # 1e8 as an arbitrary way-too-large value to capture infinities
    return(ecs)
}


setMethod("estimateDispersions", signature(cds="ExonCountSet"),
   function( cds, formula=count ~ sample + condition*exon, file=NULL, initialGuess=.01, n.cores=1, quiet=FALSE, fitDispersions=TRUE)
   {
      stopifnot(is(cds, "ExonCountSet"))
      if( all( is.na(sizeFactors(cds)) ) ){
         stop( "Estimate size factors first." )	
      }
      cds@formulas[["formula.dispersion"]] <- deparse(formula)
      # Which exons are actually testable? 
       # An all-zero gene is not testable   
      fData(cds)$testable <- rowSums(counts(cds)) > 0
      # If a gene contains less then two non-zero exons, all its exons non-testable
      fData(cds)$testable <- unlist( tapply( fData(cds)$testable, geneIDs(cds), function(x) 
         if( sum(x) > 1 ) x else rep( FALSE, length(x) ) ) ) 
      # take away genes with just one exon
      generle <- rle( as.character( geneIDs(cds) ) )
      fData(cds)$testable[which(geneIDs(cds) %in% generle$values[which( generle$lengths ==1 )])] <- FALSE

      testableGenes <- names( which( tapply( fData(cds)$testable, geneIDs(cds), any ) ) )
      if( !quiet & n.cores==1 ) {
         cat( "Dispersion estimation. (Progress report: one dot per 100 genes)\n" )
      }
      i <- 0
      if(n.cores > 1){
         if(!quiet){
            cat(sprintf("Estimating Cox-Reid exon dispersion estimates using %d cores\n", n.cores)) }
         toapply <- function(x){estimateDispersions(x, formula=formula, file=file, initialGuess=initialGuess, n.cores=1, quiet=quiet, fitDispersions=FALSE)}
         cds <- divideWork(cds, funtoapply=toapply, fattr="dispersion_CR_est", mc.cores=n.cores, testableGenes)
         cds <- fitDispersions(cds)
         cds
      }else{
         modelFrames <- sapply( testableGenes, function(gs) {
            mf <- modelFrameForGene(cds, gs)
            mf$offset <- log(mf$sizeFactor)
#            if( !quiet ) {
#               i <<- i + 1
#               if( i %% 100 == 0 )
#                  cat( "." )
#            }
            mf }, simplify=FALSE )
         
#         if( !quiet ){
#            cat( "\n  Setting up model matrices." )}
         modelmatrices <- sapply( modelFrames, function(mf){
            model.matrix( formula, data = mf )}, simplify=FALSE )

#         if( !quiet )
#            cat( "\n  Calculating initial fits." )

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
   
         fData(cds)$dispersion_CR_est <- NA_real_
   
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
   
            rows <- as.character(geneIDs(cds)) %in% genename
#        stopifnot( identical( names(disps), unname(exonIDs(cds)[rows])) ) for the momment
            fData(cds)$dispersion_CR_est[rows] <- disps 
      
            if( !quiet ) {
               i <- i + 1
               if( i %% 100 == 0 )
                  cat( "." ) }      
         }
   
#         if( !quiet & fitDispersions ){
#            cat( "\n  Fitting mean-dispersion relation" )}
   
         if( fitDispersions ){
            cds <- fitDispersions( cds )  }
   
#         if( !quiet ){
#            cat( "Finished with dispersion estimation." )}

         cds
      }
   }
)

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
      formula0 <- count ~ sample + exon + condition}
   if(is.null(formula1)){
      formula1 <- count ~ sample + exon + condition * I(exon == exonID)}
   

   # This makes sure that the formula see the 'exonID' variable used
   # below even if it the formula was supplied as parameter
   environment(formula0) <- environment()
   environment(formula1) <- environment()
   
   fam <- negative.binomial( 1 / mf$dispersion )
   fam$family = "Negative Binomial(varying)"
   
   fit0 <- try( 
         glm( formula0, data=mf, family = fam, offset = log(mf$sizeFactor), control=glm.control ),
         silent=TRUE )
   if( inherits( fit0, "try-error" ) ) {
      warning( sprintf( "Error in fit0 for gene %s: %s", gene, fit0 ) )
      return(ans) }

   for( exonID in unique(as.character(mf$exon)) ) {
      if( !nonzero[exonID] ) { 
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

testForDEU <- function( ecs, formula0=NULL, formula1=NULL, padjust=TRUE, n.cores=1, quiet=FALSE )
{
   stopifnot(is(ecs, "ExonCountSet"))
   if(sum(is.na(featureData(ecs)$dispersion))==nrow(counts(ecs))){
      stop("No dispersion parameters found, first call function estimateDispersions...\n")}

   testablegenes <- names( which( tapply( fData(ecs)$testable, geneIDs(ecs), any ) ) )
   if(n.cores > 1){
      if(!quiet){
         cat(sprintf("Testing for differential exon usage using %d cores\n", n.cores)) }
      toapply <- function(x){testForDEU(x, formula0=formula0, formula1=formula1, quiet=quiet, padjust=FALSE, n.cores=1)}
      ecs <- divideWork(ecs, funtoapply=toapply, fattr="pvalue", mc.cores=n.cores, testablegenes)
      fData(ecs)$padjust <- p.adjust(fData(ecs)$pvalue, method="BH")
    }else{
      if( !quiet )
         cat( "Testing for differential exon usage\n" )
      i <- 0
      for( genename in testablegenes ) {
         i <- i + 1
         if( !quiet & i %% 100 == 0 )
            cat( "." )
         rows <- as.character(geneIDs(ecs)) %in% genename
         res <- testGeneForDEU( ecs, genename, formula0, formula1 )
         fData(ecs)$pvalue[rows] <- res$pvalue }
      if(padjust){
         fData(ecs)$padjust <- p.adjust(fData(ecs)$pvalue, method="BH")
      }
   }
   
   ######### STORE FORMULAS IN THE OBJECTs
   if(is.null(formula0)){
      formula0 <- count~sample+exon+condition
   }
   if(is.null(formula1)){
      formula1 <- count~sample+exon+condition*I(exon==exonID)
   }
   ecs@formulas[["formula0"]] <- deparse(formula0)
   ecs@formulas[["formula1"]] <- deparse(formula1)
#   cat( "\n" )
   ecs
}

estimatelog2FoldChanges <- function(ecs, fitExpToVar="condition", n.cores=1, quiet=FALSE)
{
   stopifnot(is(ecs, "ExonCountSet"))
   if(any(is.na(sizeFactors(ecs)))){
      stop("Please estimate sizeFactors first\n")}
   if(!fitExpToVar %in% ecs@designColumns){
      stop("fitExpToVar parameter is not in the design columns, double check ecs@designColumns")}
   if(sum(is.na(featureData(ecs)$dispersion))==nrow(counts(ecs))){
      stop("No dispersion parameters found, first call function estimateDispersions...\n")}

   nms <- sort(levels(design(ecs, drop=FALSE)[[fitExpToVar]]))
   colfoldnames <- sapply(2:length(nms), function(x){
      paste("log2fold(", nms[x], "/", nms[1], ")", sep="")
   })
   
#   if(all(colfoldnames %in% colnames(fData(ecs)))){
#      stop("fData(ecs) already contains this fold changes")}

   testablegenes <- names( which( tapply( fData(ecs)$testable, geneIDs(ecs), any ) ) )
   logfold <- data.frame(matrix(ncol=length(colfoldnames), nrow=nrow(fData(ecs))))
   colnames(logfold) <- colfoldnames
   rownames(logfold) <- rownames(fData(ecs))

   if(!any(colfoldnames %in% colnames(fData(ecs)))){
      fData(ecs) <- cbind(fData(ecs), logfold)
   }

   if(n.cores==1){
      if( !quiet )
         cat( "Calculating fold changes\n" )

      frm <- as.formula(paste("count ~", fitExpToVar,  "* exon"))
      colstosteal <- paste(fitExpToVar, ":exon", sep="")
      colstopass <- c(2:length(nms))
      i <- 0
      for(gene in testablegenes){
         i <- i + 1
         if( !quiet & i %% 100 == 0 ){
            cat( "." )}
         vals <- log2(exp(t(fitAndArrangeCoefs( ecs, gene, frm=frm)[[colstosteal]])))[,colstopass]
         logfold[as.character(geneIDs(ecs)) %in% gene, colfoldnames] <- vals
      }
      fData(ecs)[,colfoldnames] <- logfold
      ecs
   }else{
      if(!quiet){
         cat(sprintf("Calculating fold changes using %d cores\n", n.cores)) }
      toapply <- function(x){estimatelog2FoldChanges(x, fitExpToVar=fitExpToVar, n.cores=1, quiet=quiet)}
      ecs <- divideWork(ecs, toapply, fattr=colfoldnames, mc.cores=n.cores, testablegenes)
      ecs
   }
}


divideWork <- function(ecs, funtoapply, fattr, mc.cores, testablegenes)
{
   if(!is.loaded("mc_fork", PACKAGE="multicore")){
      stop("Please load first multicore package or set parameter n.cores to 1...")}
   forsubset <- sort(rep(1:mc.cores, length.out=length(testablegenes)))
   subgenes <- split(testablegenes, forsubset)
   allecs <- sapply(subgenes, function(x){
      subsetByGenes(ecs, x)})
   mc.lapply <- get("mclapply", envir=getNamespace("multicore"))
   allecs <- mc.lapply(allecs, FUN=funtoapply, mc.cores=mc.cores)
   if(length(fattr) > 1){
      fData(ecs)[as.character(geneIDs(ecs)) %in% testablegenes,][,fattr] <- do.call(rbind, lapply(allecs, function(x){fData(x)[,fattr]}))
   }else{
      fData(ecs)[,fattr][as.character(geneIDs(ecs)) %in% testablegenes] <- do.call(c, sapply(allecs, function(x){fData(x)[,fattr]}))
   }
   ecs
}
