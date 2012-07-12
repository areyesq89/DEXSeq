estimateDispersionsTRT <- function( ecs, minCount=10, nCores=1 ) {
   stopifnot( inherits( ecs, "ExonCountSet" ) )
   stopifnot( "condition" %in% colnames( pData(ecs) ) )
   stopifnot( length( levels( pData(ecs)$condition ) ) == 2 )
   
   treated <- as.integer( pData(ecs)$condition ) == 2
   modelMatrix <- getModelMatrixTRT( rownames(pData(ecs)), treated )

   testable <- rowSums(counts(ecs)) >= minCount
   for( r in split( 1:nrow(ecs), geneIDs(ecs) ) ) {
      if( sum( testable[r] ) <= 1 )
         testable[r] <- FALSE
   }
   fData(ecs)$testable <- testable
   
   if( nCores == 1 ) {

      disps <- rep( NA_real_, nrow(counts(ecs)) )
      for( i in 1:length(disps) ) 
         if( fData(ecs)$testable[i] )
            try( {
               disps[i] <- 
                  estimateOneDispersionTRT( 
                     modelMatrix, 
                     getCountVectorTRT( 
                        ecs,  
                        geneIDs(ecs)[i],
                        exonIDs(ecs)[i] ),
                     sizeFactors(ecs) ) 
                if(i %% 100 == 0 )
                   cat( "." ) 
               },
            silent=FALSE) 
      attr( disps, "TRT" ) <- TRUE
      fData(ecs)$dispBeforeSharing <- disps

   } else {

      testableGenes <- as.character(unique(fData(ecs)[which(fData(ecs)$testable),]$geneID))
      ecs <- divideWork( 
         ecs, 
         funtoapply = estimateDispersionsTRT,
         fattr = "dispBeforeSharing", 
         mc.cores = nCores, 
         testableGenes )

   }
   ecs
}

testForDEU.TRT <- function( ecs, nCores=1 ) {
   stopifnot( inherits( ecs, "ExonCountSet" ) )
   if( all(is.na(featureData(ecs)$dispersion)) )
      stop("No dispersion values found, call function fitDispersionFunction first.")
   
   sf <- rep( sizeFactors( ecs ), 2 )
   disps <- fData(ecs)$dispersion
   treated <- as.integer( pData(ecs)$condition ) == 2
   mm1 <- getModelMatrixTRT( rownames(pData(ecs)), treated )
   mm0 <- mm1[,-ncol(mm1)]
   
   if( nCores == 1 ) {

      pvals <- rep( NA_real_, nrow(counts(ecs)) )
      for( i in 1:length(pvals) ) 
         if( fData(ecs)$testable[i] ) {
            cv <- getCountVectorTRT( ecs, geneIDs(ecs)[i], exonIDs(ecs)[i] )
            try( 
               {
                  fit1 <- statmod::glmnb.fit( mm1, cv, disps[i], log(sf) )
                  fit0 <- statmod::glmnb.fit( mm0, cv, disps[i], log(sf) )
                  pvals[i] <- 1 - pchisq( deviance(fit0) - deviance(fit1), 1 )
               if(i %% 100 == 0 )
                   cat( "." ) 
               },
               silent=TRUE)
         }
      fData(ecs)$pvalue <- pvals
   
   } else{

      testableGenes <- as.character(unique(fData(ecs)[which(fData(ecs)$testable),]$geneID))
      ecs <- divideWork( 
         ecs, 
         funtoapply = testForDEU.TRT,
         fattr = "pvalue", 
         mc.cores = nCores, 
         testableGenes )

   }

   fData(ecs)$padjust <- p.adjust( fData(ecs)$pvalue, method="BH" )

   ecs
}


# internal functions

getCountVectorTRT <- function( ecs, geneID, exonID ) {
   stopifnot( inherits( ecs, "ExonCountSet" ) )
   ctfg <- countTableForGene( ecs, geneID )
   w <- which( rownames(ctfg) == exonID )   
   stopifnot( length(w) == 1 )
   c( ctfg[ w, ], colSums( ctfg[ -w, , drop=FALSE ] ) )
}

getModelMatrixTRT <- function( sampleNames, treated ) {
   mf <- data.frame(
      sample = rep( sampleNames, 2 ),
      treated = rep( ifelse( treated, "yes", "no" ), 2 ),
      what = rep( c( "this", "others" ), each = length(sampleNames) ) )
   
   mm <- model.matrix( ~ 0 + sample + what + treated:what, mf )
   mm <- mm[ , colnames(mm) != "whatothers:treatedyes" ]
   rownames(mm) <- paste( 
      rep( sampleNames, 2), 
      rep( c( "this", "other" ), each = length(sampleNames) ),
      sep = ":" )
   mm
}

estimateOneDispersionTRT <- function( modelMatrix, countVector, 
      sizeFactors, dispInitialGuess = .5 ) {
   
   fit1 <- statmod::glmnb.fit( modelMatrix, countVector, 
                      dispInitialGuess, rep( log(sizeFactors), 2 ) ) 
   
   exp( optimize( 
      function( logalpha ) 
         -profileLogLikelihood( 
            exp(logalpha), 
            modelMatrix,
            countVector, 
            fitted.values(fit1) ), 
      log( c( 1e-5, 1e3 ) ) )$minimum )
}

#profileLogLikelihood <- DEXSeq:::profileLogLikelihood
