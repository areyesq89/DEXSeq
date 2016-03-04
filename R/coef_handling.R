arrangeCoefs <- function( frm, mf, mm = model.matrix( frm, mf ), fit = NULL, insertValues = TRUE ) {

   if( any( attr( mm, "contrasts" ) != "contr.treatment" ) )
      stop( "Can only deal with standard 'treatment' contrasts." )   # Do I need this?
   if( is.null(fit) & insertValues )
      stop( "If fit==NULL, returnCoefValues must be FALSE" )
   if( !is.null(fit) )
      stopifnot( all( colnames(mm) == names(coefficients(fit)) ) )

   fctTbl <- attr( terms(frm), "factors" )

   coefIndicesList <- 
   lapply( seq_len(ncol(fctTbl)), function( fctTblCol ) {
      termName <- colnames(fctTbl)[ fctTblCol ]
      varsInTerm <- stringr::str_split( termName, stringr::fixed(":") )[[1]] 
      stopifnot( all( fctTbl[ varsInTerm, fctTblCol ] == 1 ) )
      stopifnot( sum( fctTbl[ , fctTblCol ] ) == length( varsInTerm ) )
      coefNames <- colnames(mm)[ attr( mm, "assign" ) == fctTblCol ]
      lvlTbl <- stringr::str_match( coefNames, 
         stringr::str_c( "^", stringr::str_c( varsInTerm, "([^:]*)", collapse=":" ), "$" ) )[ , -1, drop=FALSE ]
      stopifnot( ncol(lvlTbl) == length( varsInTerm ) )
      stopifnot( nrow(lvlTbl) == length( coefNames ) )
      if( !all( sapply( varsInTerm, function(v) is.factor(mf[[v]]) | is.character(mf[[v]]) ) ) )
         stop( "Non-factor in model frame" )

      varLevels <- lapply( varsInTerm, function(v) levels( factor( mf[[v]] ) ) ) 
      coefIndices <- array( NA_character_, dim = sapply( varLevels, length ), dimnames = varLevels )
      names( dimnames( coefIndices ) ) <- varsInTerm

      for( i in seq_len( nrow(lvlTbl) ) )
         coefIndices <- do.call( `[[<-`, c( quote(coefIndices), as.list( lvlTbl[ i, ] ), coefNames[i] ) )

      coefIndices
   } )
   names( coefIndicesList ) <- colnames( fctTbl )

   if( attr( terms(frm), "intercept" ) ) {
      a <- array( c( `(Intercept)` = "(Intercept)" ) )
      dimnames(a) <- list( `(Intercept)` = c( "(Intercept)" ) )
      coefIndicesList <- c( list( `(Intercept)` = a ), coefIndicesList )
   }

   if( !insertValues )
      ans <- coefIndicesList
   else
      ans <- lapply( coefIndicesList, function(coefIndices) {
         a <- ifelse( is.na(coefIndices), 0, coefficients(fit)[ coefIndices ] )
         attr( a, "variables" ) <- attr( coefIndices, "variables" )
         a } )
      
   lapply( ans, function(x) 
      if( is.array(x) ) 
         x 
      else { 
         y <- array( x, dim=length(x) )
         attr( y, "variables" ) <- attr( x, "variables" )
         dimnames(y) <- list( names(x) )
         y } )
}

apply2 <- function( X, MARGIN, FUN, ... ) {
   if( length(MARGIN) > 0 ) 
      apply( X, MARGIN, FUN, ... ) 
   else 
      FUN( X, ... ) }

balanceExons <- function( coefs, dispersions ) {
   stopifnot( any( sapply( coefs, function(x) 
      identical( names(dimnames(x)), "(Intercept)" ) ) ) )
   termsWithExon <- sapply( coefs, function(x) "exon" %in% names(dimnames(x)) )
   meanMainEffect <- sum( sapply( coefs[!termsWithExon], mean, na.rm=TRUE ) )
   meanExonEffects <- rowSums( sapply( coefs[termsWithExon], function(x) 
      apply2( x, "exon", mean, na.rm=TRUE ) ) )

   meanExonFittedValue <- exp( meanMainEffect + meanExonEffects )

   exonWeights <-  1 / ( dispersions + 1 / meanExonFittedValue )

   shifts <- lapply( coefs[termsWithExon], function(x) { 
      nonExonDims <- which(  names(dimnames(x)) != "exon" )
      list(
         vars = names(dimnames(x))[ nonExonDims ],
         wmeans = apply2( x, nonExonDims, weighted.mean, exonWeights) ) } )

   lapply( coefs, function(x) {
      nonExonVars <- names(dimnames(x))[ names(dimnames(x)) != "exon" ]
      if( identical( nonExonVars, "(Intercept)" ) )
         whichShift <- which( sapply( shifts, function(xx) length( xx$vars ) == 0 ) )
      else
         whichShift <- which( sapply( shifts, function(xx) identical( xx$vars, nonExonVars ) ) )
      if( length( whichShift ) == 0 )
         return( x )
      if( length( whichShift ) > 1 )
         stop( "Confused about selecting shift." )
      if( "exon" %in% names(dimnames(x)) )
         x - shifts[[ whichShift ]]$wmeans
      else
         x + shifts[[ whichShift ]]$wmeans
    } )
}         

fitAndArrangeCoefs <- function( frm = count ~ condition * exon, balanceExons = TRUE, mf, fitExpToVar)
{
   if( length(levels(mf$exon)) <= 1 )
      return( NULL )
   mm <- model.matrix( frm, mf )
   fit <- try(
              glmnb.fit(mm,
                        mf$count,
                        dispersion=mf$dispersion,
                        offset=log(mf$sizeFactor),
                        tol=0.1),
              silent=TRUE)
   
   if( is( fit, "try-error" ) ){
       return( NULL )
   }

   if( is( mf[[fitExpToVar]], "numeric" ) ){
       coefs <- fit$coefficients
       attributes(coefs)$fitType <- "numeric"
   }else{
       coefs <- arrangeCoefs( frm, mf, mm, fit )
       if( balanceExons ) {
           coefs <- balanceExons( coefs, tapply( mf$dispersion, mf$exon, `[`, 1 ) )
       }
       attributes(coefs)$fitType <- "factor"
   }
   coefs
}

getEffectsForPlotting <- function( coefs, groupingVar = "condition", averageOutExpression=FALSE, frm, mf )
{
    if( attributes(coefs)$fitType == "factor" ){
        groupingExonInteraction <- which( sapply( coefs, function(x) 
            all( c( groupingVar, "exon") %in% names(dimnames(x)) ) & length(dim(x)) == 2 ) ) 
        fittedValues <- coefs[[ groupingExonInteraction ]]
        if( names(dimnames(fittedValues))[1] == "exon" )
            fittedValues <- t( fittedValues )
        stopifnot( identical( names(dimnames(fittedValues)), c( groupingVar, "exon" ) ) )
        for( x in coefs[ -groupingExonInteraction ] ) {
            if( all( c( groupingVar, "exon") %in% names(dimnames(x)) ) )
                stop( "Cannot yet deal with third-order terms." )
            if( !any( c( groupingVar, "exon") %in% names(dimnames(x)) ) ) {
                fittedValues <- fittedValues + mean( x )
            } else if( averageOutExpression & identical( names(dimnames(x)), groupingVar ) ) {
                fittedValues <- fittedValues + mean( x )
            } else if( groupingVar %in% names(dimnames(x)) ) {
                groupMeans <- apply2( x, groupingVar, mean )
                stopifnot( identical( names(groupMeans), dimnames(fittedValues)[[1]] ) )
                fittedValues <- fittedValues + groupMeans
            } else if( "exon" %in% names(dimnames(x)) ) {
                exonMeans <- apply2( x, "exon", mean )
                fittedValues <- t( t(fittedValues) + exonMeans )
            } else {
                print( x )
                stop( "Unexpected term encountered." )
            }
        }
        return( fittedValues )
    }else{
        stopifnot( "(Intercept)" %in% names(coefs) )
        stopifnot( "exonthis" %in% names(coefs) )
        allVars <- all.vars(frm)
        continuousVar <- allVars[!allVars %in% c("count", "exon")]
        interactionCoefName <- paste0( continuousVar, ":exonthis" )
        stopifnot( interactionCoefName %in% names(coefs) )
        mf[[continuousVar]]
        predictors <- unique( mf[[continuousVar]] )
        if( averageOutExpression ){
            fittedValues <- coefs["exonthis"] + coefs[interactionCoefName]*predictors
        }else{
            fittedValues <- coefs["(Intercept)"] + coefs["exonthis"] +
                (coefs[continuousVar] + coefs[interactionCoefName])*predictors
        }
        fittedValues <- matrix(fittedValues, ncol=1)
        colnames(fittedValues) <- "this"
#        rownames(fittedValues) <- sprintf("%s=%s", continuousVar, predictors)
        rownames(fittedValues) <- as.character(predictors)
        return( fittedValues )
    }
}


modelFrameSM <- function(object)
{
    mfSmall <- as.data.frame( colData(object) )
    mfSmall$exon <- relevel(mfSmall$exon, "others")
    mfSmall$dispersion <- NA
    mfSmall$count <- NA
    mfSmall
}

getEffectsForGeneBM <- function(geneID, groups, notNAs, countsAll,
                                disps, features, mf, frm, numsamples,
                                fitExpToVar, averageOutExpression=TRUE)
{
    rt <- groups %in% geneID & notNAs
    if( sum(rt) < 2 ){ return(NULL) }
    countsThis <- countsAll[rt,]
    rownames(countsThis) <- gsub("\\S+:", "", rownames(countsThis))
    dispsThis <- disps[rt]
    names(dispsThis) <- features[rt]
    numexons <- sum(rt)
    newMf <- mf[as.vector( sapply( split( seq_len(nrow(mf)), mf$sample ), "[", seq_len( numexons ) ) ),]
    newMf$exon <- factor( rep( features[rt], numsamples ) )
    for (i in seq_len(nrow(newMf))) {
       newMf[i, "dispersion"] <- dispsThis[as.character(newMf[i, "exon"])]
       newMf[i, "count"] <- countsThis[as.character(newMf[i, "exon"]), as.character(newMf[i, "sample"])]
    }
    newMf <- droplevels(newMf)
    coefficients <- fitAndArrangeCoefs( frm, balanceExons = TRUE, mf=newMf, fitExpToVar=fitExpToVar)
    if (is.null(coefficients)){
       return(coefficients)
    }
    ret <- t( getEffectsForPlotting(coefficients, averageOutExpression = averageOutExpression, 
        groupingVar = fitExpToVar, frm=frm, mf=newMf))
    rownames(ret) <- paste(geneID, rownames(ret), sep = ":")
    return(ret)
}

getEffectsForExonsSM <- function(index, frm, countsAll, disps,
                                 mfSmall, averageOutExpression=TRUE,
                                 fitExpToVar)
{
    mfSmall$count <- countsAll[index,]
    mfSmall$dispersion <- disps[index]
    getEffectsForPlotting(
        fitAndArrangeCoefs( frm, mf=mfSmall, balanceExons=FALSE, fitExpToVar=fitExpToVar),
            averageOutExpression=averageOutExpression,
            groupingVar=fitExpToVar, frm, mfSmall )[,"this"]
}

getEffectsForGene <- function( geneID, object, maxRowsMF, fitExpToVar)
{
    rt <- object$groupID %in% geneID
    sampleData <- object@sampleData
    if( is(object@sampleData[[fitExpToVar]], "numeric") ){
        maxRowsMF <- 0
    }
    numsamples <- nrow(object@sampleData)
    numexons <- sum(rt)
    featuresInGene <- object$featureID[rt]
    dispersions <- object$dispersion[rt]
    dispersions[is.na(dispersions)] <- 1e-08
    frm <- as.formula(paste("count ~", fitExpToVar, "* exon"))
    bigFlag <- numsamples*numexons < maxRowsMF
    if( bigFlag ){
        mf <- object@modelFrameBM
        mf <- mf[as.vector(sapply(split(seq_len(nrow(mf)), mf$sample), 
            "[", seq_len(numexons))), ]
        mf$exon <- factor(rep(featuresInGene, nrow(sampleData)))
        counts <- object$countData[rt,]
        rownames(counts) <- gsub("\\S+:", "", rownames(counts))
        names(dispersions) <- object$featureID[rt]
        for (i in seq_len(nrow(mf))) {
            mf[i, "dispersion"] <-
                dispersions[as.character(mf[i, "exon"])]
            mf[i, "count"] <-
                counts[as.character(mf[i, "exon"]), as.character(mf[i, "sample"])]
        }
        mf <- droplevels(mf)
        coefs <- fitAndArrangeCoefs(frm, balanceExons=TRUE, mf=mf, fitExpToVar=fitExpToVar)
        if( is.null(coefs ) ){
            return()
        }
        splicing <- t(getEffectsForPlotting( coefs, groupingVar=fitExpToVar, averageOutExpression=TRUE, frm=frm, mf=mf))
        expression <- t(getEffectsForPlotting( coefs, groupingVar=fitExpToVar, averageOutExpression=FALSE, frm=frm, mf=mf))
        rownames(splicing) <- sprintf("%s:%s", geneID, rownames(splicing))
        rownames(expression) <- rownames(splicing)
        list( expression=expression, splicing=splicing )
    }else{
        mf <- object@sampleData
        mf <- rbind( data.frame(mf, exon="this"), data.frame(mf, exon="others"))
        mf$exon <- relevel( mf$exon, "others" )
        countsThis <- object$countData[rt,]
        countsOthers <- sapply( rownames( countsThis ),
                               function(x){
                                   colSums(countsThis[!rownames(countsThis) %in% x,,drop=FALSE])
                               })
        countsOthers <- t(countsOthers)
        stopifnot(all(rownames(countsThis) ==  rownames(countsOthers)))
        effects <- lapply( seq_len(numexons), function(x){
                   mf$count <- c( countsThis[x,], countsOthers[x,])
                   mf$dispersion <- dispersions[x]
                   coefs <- fitAndArrangeCoefs(frm, balanceExons=FALSE, mf=mf, fitExpToVar=fitExpToVar)
                   if( is.null(coefs) ){
                       return(NULL)
                   }
                   splicing <- getEffectsForPlotting( coefs, groupingVar=fitExpToVar, averageOutExpression=TRUE, frm=frm, mf=mf)[,"this"]
                   expression <- getEffectsForPlotting( coefs, groupingVar=fitExpToVar, averageOutExpression=FALSE, frm=frm, mf=mf)[,"this"]
                   list(splicing=splicing, expression=expression)
               })
        names(effects) <- rownames(object)[rt]
        splicing <- t(sapply(effects, "[[", "splicing"))
        expression <- t( sapply(effects, "[[", "expression" ))
        list( expression=expression, splicing=splicing )
    }
}
