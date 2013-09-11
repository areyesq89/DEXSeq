setMethod("estimateSizeFactors", signature(object="ExonCountSet"),
   function( object ){
      cds <- object
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


profileLogLikelihood <- function( disp, mm, y, muhat )
# FIXME: Change function name; this is a conditional likelihood, not
# a profile likelihood. 
{
   # calculate the log likelihood:
   if(length(disp) != length(y)){
      disp <- rep(disp, length(y))
   }

   ll <- sum( sapply( seq(along=y), function(i)
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
   cr <- sum( log( abs( diag( qrres$qr )[ seq_len(qrres$rank) ] ) ) )

   # return the profile log likelihood:
   ll - cr
}


fitDispersionFunction <- function( ecs )
{
   stopifnot(is(ecs, "ExonCountSet"))
   if(all(is.na(fData(ecs)$dispBeforeSharing))){
      stop("no CR dispersion estimations found, please first call estimateDispersions function")
   }
   means <- colMeans( t(counts(ecs))/sizeFactors(ecs) )
   disps <- fData(ecs)$dispBeforeSharing
   coefs <- c( .1, 1 )
   iter <- 0
   while(TRUE) {
      residuals <- disps / ( coefs[1] + coefs[2] / means )
      good <- which((residuals > 1e-4) & (residuals < 15))
      mm <- model.matrix(disps[good] ~ I(1/means[good]))
      fit <- try(glmgam.fit(mm, disps[good], coef.start=coefs), silent=TRUE)
      if(inherits(fit, "try-error")){
         stop("Failed to fit the dispersion function\n")
      }
      oldcoefs <- coefs
      coefs <- coefficients(fit)
      if(coefs[1] < 0){
         coefs[1] <- 0
         warning("Negative intercept value in the dispersion function, it will be set to 0. Check fit diagnostics plot section from the vignette.")
         break
      }
      if( sum( log( coefs / oldcoefs )^2 ) < .005 )
         break
      iter <- iter + 1
      if( iter > 10 ) {
         warning( "Dispersion fit did not converge." )
         break }
    }
    ecs@dispFitCoefs <- coefs
    fData(ecs)$dispFitted <- ecs@dispFitCoefs[1] + ecs@dispFitCoefs[2] / colMeans( t(counts(ecs))/sizeFactors(ecs) )
    fData(ecs)$dispersion <- pmin(
       pmax(
          fData(ecs)$dispBeforeSharing,
          fData(ecs)$dispFitted,
          na.rm = TRUE ),
          1e8 )   # 1e8 as an arbitrary way-too-large value to capture infinities
    return(ecs)
}





estimatelog2FoldChanges <- function(ecs, fitExpToVar="condition", denominator="", getOnlyEffects=FALSE, averageOutExpression=TRUE, nCores=1, quiet=FALSE, file="")
{
   stopifnot(is(ecs, "ExonCountSet"))
   if(any(is.na(sizeFactors(ecs)))){
      stop("Please estimate sizeFactors first\n")}
   if(!fitExpToVar %in% ecs@designColumns){
      stop("fitExpToVar parameter is not in the design columns, double check ecs@designColumns")}
   if(sum(is.na(featureData(ecs)$dispersion))==nrow(counts(ecs))){
      stop("No dispersion parameters found, first call function estimateDispersions...\n")}

   frm <- as.formula(paste("count ~", fitExpToVar,  "* exon"))
   testablegenes <- as.character(unique(fData(ecs)[which(fData(ecs)$testable),]$geneID))

   geteffects <- function(geneID){
     coefficients <- fitAndArrangeCoefs(ecs, geneID=geneID, frm, balanceExons=TRUE)
     if( is.null( coefficients ) ){
        return(coefficients)
     }
     ret <- t(getEffectsForPlotting(coefficients, averageOutExpression=averageOutExpression, groupingVar=fitExpToVar))
     rownames(ret) <- paste(geneID, rownames(ret), sep=":")
     return(ret)
   }

   if( nCores > 1 ){
      if(!is.loaded("mc_fork", PACKAGE="parallel")){
      stop("Please load first parallel package or set parameter nCores to 1...")}
      alleffects <- parallel:::mclapply( testablegenes, function(x){ geteffects(x) }, mc.cores=nCores )
    }else{
      alleffects <- lapply( testablegenes, function(x){geteffects(x)})
    }

    names(alleffects) <- testablegenes
    alleffects <- do.call(rbind, alleffects)
    alleffects <- vst(exp( alleffects ), ecs)
    toadd <- matrix(NA, nrow=nrow(ecs), ncol=ncol(alleffects))
    rownames(toadd) <- featureNames(ecs)

    if( getOnlyEffects ){
       colnames(toadd) <- colnames(alleffects)
       toadd[rownames(alleffects), colnames(alleffects)] <- alleffects
     }else{
       if( denominator == "" ){
          denominator <- as.character(design(ecs, drop=FALSE)[[fitExpToVar]][1])
       }
       stopifnot( any( colnames(alleffects) %in% denominator ) )
       denoCol <- which(colnames(alleffects) == denominator)
       alleffects <- log2(alleffects / alleffects[,denoCol])
       colnames(alleffects) <- sprintf("log2fold(%s/%s)", colnames(alleffects), denominator)
       colnames(toadd) <- colnames(alleffects)
       alleffects <- alleffects[,-denoCol, drop=FALSE]
       toadd <- toadd[,-denoCol, drop=FALSE]
       toadd[rownames(alleffects), colnames(alleffects)] <- alleffects
     }

    fData(ecs) <- cbind(fData(ecs), toadd)
    ecs
}


divideWork <- function(ecs, funtoapply, fattr, mc.cores, testablegenes)
{
#   if(!suppressMessages(suppressWarnings(require("parallel")))){
 #     stop("parallel package not found...")}
   if(!is.loaded("mc_fork", PACKAGE="parallel")){
     stop("Please load first parallel package or set parameter nCores to 1...")}
   stopifnot(mc.cores>=1)

   subgenes <- split(testablegenes, seq(along=testablegenes) %% mc.cores)
   allecs <- lapply(subgenes, function(x) subsetByGenes(ecs, x) )
   allecs <- parallel::mclapply(allecs, FUN=funtoapply, mc.cores=mc.cores)

   for(j in seq(along=allecs)) {
      rownam <-  rownames(fData(allecs[[j]]))
      fData(ecs)[ rownam, fattr] <-  fData(allecs[[j]])[, fattr]
   }

   ecs
}

modelFrameForGene <- function( ecs, geneID, onlyTestable=FALSE) {
   stopifnot( is( ecs, "ExonCountSet") )
   if( onlyTestable & any(colnames(fData(ecs)) %in% "testable")){
      rows <- geneIDs(ecs) == geneID & fData(ecs)$testable
   }else{
      rows <- geneIDs(ecs) == geneID
   }

   numExons <- sum(rows)
   exonCol <- rep( factor( exonIDs(ecs)[rows], levels=exonIDs(ecs)[rows] ), ncol( counts(ecs) ) )
   modelFrame <- data.frame(
      sample = rep( factor( colnames( counts(ecs) ) ), each = numExons ),
      exon = exonCol,
      sizeFactor = rep( sizeFactors(ecs), each = numExons ) )
   for( cn in colnames( design(ecs,drop=FALSE) ) )
      modelFrame[[cn]] <- factor(rep( design(ecs,drop=FALSE)[[cn]], each=numExons ), levels=sort(levels(design(ecs,drop=FALSE)[[cn]] )))
   modelFrame$dispersion <- fData(ecs)$dispersion[ rows ][
      match( modelFrame$exon, exonIDs(ecs)[rows] ) ]
   modelFrame$count <- as.vector( counts(ecs)[rows,] )
   attr( modelFrame, "geneID" ) <- geneID
   modelFrame
}
