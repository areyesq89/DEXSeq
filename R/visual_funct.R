plotDEXSeq <- function(ecs, geneID, FDR=0.1, coefficients=TRUE, norCounts=FALSE, expression=TRUE, displayTranscripts=FALSE, names=FALSE, legend=FALSE, color=NULL, ...){
	stopifnot(is(ecs, "ExonCountSet"))
	if(any(is.na(sizeFactors(ecs)))){
		stop("Please estimate sizeFactors first\n")
	}

	if(coefficients==FALSE & norCounts==FALSE){
		stop("Please indicate what would you like to plot")
	}
	
	rt<-which(featureData(ecs)$geneID==geneID)
	rango <- 1:length(rt)
	count <- t(t(counts(ecs)[rt,])/sizeFactors(ecs))

	if(sum(count) == 0){
		warning("No read counts falling in this gene, there is nothing to plot.")
		return()
	}
	vst <- function(x)
      ( 2 / ( sqrt(ecs@dispFitCoefs[1]) ) ) * 
         log( 2 * ecs@dispFitCoefs[1] * sqrt(x) + 
              2 * sqrt( ecs@dispFitCoefs[1] * ( ecs@dispFitCoefs[2] + 1 + ecs@dispFitCoefs[1] * x ) ) ) -
      ( 2 / ( sqrt(ecs@dispFitCoefs[1]) ) ) * 
         log( 2 * sqrt( ecs@dispFitCoefs[1] * ( ecs@dispFitCoefs[2] + 1 ) ) )
	
	numcond<-length(unique(design(ecs, drop=FALSE)$condition))
	numexons<-nrow(count)
		
	##### prepare plotting values coefficients #######
	if(FDR>1|FDR<0){stop("FDR has to be a numeric value between 0 - 1")}
	
	each <- featureData(ecs)$padjust[as.character(featureData(ecs)$geneID) %in% geneID]

	##### reescale values from positions in the genome to 0 to 1 for plotting
	if(!any(is.na(featureData(ecs)$start))){
		sub<-data.frame(start=featureData(ecs)$start[rt], end=featureData(ecs)$end[rt], chr=featureData(ecs)$chr[rt], strand=featureData(ecs)$strand[rt])
		rel<-(data.frame(sub$start, sub$end))-min(sub$start)	
		rel<-rel/max(rel[,2])
	##### determine the layout #######  depending on the information available-- values tested empirically
		if(displayTranscripts==TRUE & !is.null(featureData(ecs)$transcripts)){
			transcripts <- sapply(featureData(ecs)$transcripts[featureData(ecs)$geneID %in% geneID], function(x){strsplit(x, ";")})
			trans <- Reduce(union, transcripts)
			mat <- 1:(3+min(length(trans), 45)) ## max support from transcripts is 45, which seems to be the max for the layout supported by graphics
			hei<-c(8, 1, 1.5, rep(1.5, min(length(trans), 45)))
		}else{
			mat<-1:3
			hei<-c(5, 1, 1.5)
		}
	
		if(sum(c(coefficients, norCounts))>1){
			hei <- c(hei[1], hei)
			mat <- c(mat, length(mat)+1)	
		}
		hei <- c(hei, .2)
		mat <- c(mat, length(mat)+1)
		layout(matrix(mat), height=hei)
		par(mar=c(2, 4, 4, 2))
	}else if(coefficients & norCounts){
		par(mfrow=c(2,1))
	}
	##### determine colors if not provided by user
	if(is.null(color)){
		color<-rgb(colorRamp(c("#D7191C", "#FFFFBF", "#2B83BA"))(seq(0, 1, length.out=numcond)), max=255, alpha=175)
	}

	names(color) <- sort(levels(design(ecs, drop=FALSE)$condition))

	if(coefficients){
		plot.new()
		if(sum(is.na(featureData(ecs)$dispersion))==nrow(counts(ecs))){
			stop("No dispersion parameters found, first call function estimateDispersions...\n")		
		}
		if( expression ){
			mtext("Fitted expression", side=2, adj=0.5, padj=1, line=1.5, outer=FALSE, ...)
		}else{
			mtext("Fitted splicing", side=2, adj=0.5, padj=1, line=1.5, outer=FALSE, ...)
		}

		coeff <- as.data.frame( t( getEffectsForPlotting( fitAndArrangeCoefs( ecs, geneID ), averageOutExpression=!expression) ) )
		coeff <- exp(coeff)

		ylimn <- c(min(coeff, na.rm=TRUE), max(coeff, na.rm=TRUE))
		coeff <- vst( coeff  )
		plot.window(xlim=c(0, 1), ylim=c(min(coeff), max(coeff)))
		minlog10 <- floor( log10( 1/ncol(counts(ecs)) ) )
		maxlog10 <- ceiling( log10( ylimn[2] ) )
		ticks <- 10^seq( minlog10, maxlog10 )
		decade_lengths <- ( vst(ticks)[ 2 : length(ticks) ] - vst(ticks)[ 1 : (length(ticks)-1) ] ) /
		   ( vst( ylimn[2] ) - vst( ylimn[1] ) )
		ticks <- c( 0, ticks[ min( which( decade_lengths > .1 ) ) : length(ticks) ] )
		axis( 2, at=vst(ticks), labels=ticks, las=2, pos=0, ...)

		for( i in minlog10 : (maxlog10-1) ) {
		   decade_length <- ( vst( 10^(i+1) ) - vst( 10^i ) ) / ( vst( ylimn[2] ) - vst( ylimn[1] ) )
		   if( decade_length > .1 ) {
		      axis( 2, at = vst( 1:9 * 10^i ), labels = FALSE, las=2, tcl=-.25, pos=0, ...)
		   }
		   if( decade_length > .3 & decade_length <= .6) {
		      axis( 2, at = vst( c(2,3,5) * 10^i ), labels = ( c(2,3,5) * 10^i ), las=2, tcl=-.25, pos=0, ...)
		   } else if( decade_length > .6 ) {
		      axis( 2, at = vst( c(1.5,2:9) * 10^i ), labels = ( c(1.5,2:9) * 10^i ), las=2, tcl=-.25, pos=0, ...)
		   }
		}
		axis( 2, at=vst(0), labels=0, las=2, pos=0, ...)

		intervals<-(0:nrow(coeff))/nrow(coeff)
		middle <- apply(cbind(intervals[rango], 
			(intervals[rango+1]-((intervals[rango+1])-intervals[rango])*0.2)), 
			1, median)
		coeff <- rbind(coeff, NA)
		for(j in 1:ncol(coeff)){
			segments(intervals[rango], 
			         coeff[rango,j], intervals[rango+1]-((intervals[rango+1]-intervals[rango])*0.2), 
 				 coeff[rango,j], 
				 col=color[sort(levels(design(ecs, drop=FALSE)$condition))[j]], ...)  #### line with the y level
			segments(intervals[rango+1]-((intervals[rango+1]-intervals[rango])*0.2), 
				 coeff[rango,j], 
				 intervals[rango+1], 
				 coeff[rango+1,j], 
				 col=color[sort(levels(design(ecs, drop=FALSE)$condition))[j]], lty="dotted", ...)  #### line joining the
		}
		abline(v=middle[rango], lty="dotted", col="lightgrey")
		axis(1, at=middle[1:length(rt)], labels=featureData(ecs)$exonID[rt], ...)
	}

	if(norCounts){
		plot.new()
		count <- countTableForGene(ecs, geneID, normalized=TRUE)
		ylimn <- c(min(count, na.rm=TRUE), max(count, na.rm=TRUE))
		count <- vst( count )
		plot.window(xlim=c(0, 1), ylim=c(min(count), max(count)))
		minlog10 <- floor( log10( 1/ncol(counts(ecs)) ) )
		maxlog10 <- ceiling( log10( ylimn[2] ) )
		ticks <- 10^seq( minlog10, maxlog10 )
		decade_lengths <- ( vst(ticks)[ 2 : length(ticks) ] - vst(ticks)[ 1 : (length(ticks)-1) ] ) /
		   ( vst( ylimn[2] ) - vst( ylimn[1] ) )
		ticks <- c( 0, ticks[ min( which( decade_lengths > .1 ) ) : length(ticks) ] )
		axis( 2, at=vst(ticks), labels=ticks, las=2, pos=0, ...)

		for( i in minlog10 : (maxlog10-1) ) {
		   decade_length <- ( vst( 10^(i+1) ) - vst( 10^i ) ) / ( vst( ylimn[2] ) - vst( ylimn[1] ) )
		   if( decade_length > .1 ) {
		      axis( 2, at = vst( 1:9 * 10^i ), labels = FALSE, las=2, tcl=-.25, pos=0, ...)
		   }
		   if( decade_length > .3 & decade_length <= .6) {
		      axis( 2, at = vst( c(2,3,5) * 10^i ), labels = ( c(2,3,5) * 10^i ), las=2, tcl=-.25, pos=0, ...)
		   } else if( decade_length > .6 ) {
		      axis( 2, at = vst( c(1.5,2:9) * 10^i ), labels = ( c(1.5,2:9) * 10^i ), las=2, tcl=-.25, pos=0, ...)
		   }
		}
		axis( 2, at=vst(0), labels=0, las=2, pos=0, ...)

		intervals<-(0:nrow(count))/nrow(count)
		middle <- apply(cbind(intervals[rango], (intervals[rango+1]-((intervals[rango+1])-intervals[rango])*0.2)), 1, median)
		count<-rbind(count, NA)

		for(j in 1:ncol(count)){
			segments(intervals[rango], 
			         count[rango,j], 
			         intervals[rango+1]-((intervals[rango+1]-intervals[rango])*0.2), 
			         count[rango,j], 
			         col=color[as.character(design(ecs, drop=FALSE)$condition)[j]], ...)  #### line with the y level
			segments(intervals[rango+1]-((intervals[rango+1]-intervals[rango])*0.2), 
			         count[rango,j], 
				 intervals[rango+1], 
				 count[rango+1,j], 
				 col=color[as.character(design(ecs, drop=FALSE)$condition)[j]], lty="dotted", ...)  #### line joining the y levels
		}

		abline(v=middle[rango], lty="dotted", col="lightgrey")
		mtext("Normalized counts", side=2, adj=0.5, padj=1, line=1.5, outer=FALSE, ...)
		axis(1, at=middle[1:length(rt)], labels=featureData(ecs)$exonID[rt], ...)

	}

	########### plot the gene model ########## just if transcript information available
	if( !any(is.na(featureData(ecs)$start)) ){
		par(mar=c(0, 4, 0, 2))
		plot.new()
		exoncol<-ifelse(each<=FDR, "#8B0000", "dark green")
		exoncol[is.na(exoncol)]<-"black"
		segments(apply((rbind(rel[rango,2], rel[rango, 1])), 2, median), 0, apply(rbind(intervals[rango], intervals[rango+1]-((intervals[rango+1]-intervals[rango])*0.2)), 2, median), 1, col=exoncol)
		par(mar=c(1.5, 4, 0, 2))
		plot.new()
		plot.window(xlim=c(min(sub$start), max(sub$end)), ylim=c(0, 1))
		rect(sub[rango,1], 0, sub[rango,2], 1, col=exoncol)
		segments(sub[rango,2], 0.5, apply(rbind(sub[rango, 2], sub[rango+1, 1]), 2, median), 0.65)
		segments(apply(rbind(sub[rango, 2], sub[rango+1, 1]), 2, median), 0.65, sub[rango+1,1], 0.5)
		if(names){
		mtext("Gene Model", side=2, adj=0.5, padj=1, line=0, outer=FALSE, cex=1, las=1)
		}
	
		if(!is.null(featureData(ecs)$transcripts)){
		##### plot the transcripts #######
			if(displayTranscripts){
			for(i in 1:min(length(trans), 45)){
			    logicexons <- sapply(transcripts, function(x){length(which(x==trans[i]))})
			    tr<-data.frame(featureData(ecs)$start[rt][logicexons==1], featureData(ecs)$end[rt][logicexons==1])
			    plot.new()
		            plot.window(xlim=c(min(sub$start), max(sub$end)), ylim=c(0, 1))
        			rect(tr[rango,1], 0, tr[rango,2], 1)
				segments(tr[rango,2], 0.5, apply(rbind(tr[rango,2], tr[rango+1, 1]), 2, median), 0.65)
				segments(apply(rbind(tr[rango,2], tr[rango+1, 1]), 2, median), 0.65, tr[rango+1,1], 0.5)
			    if(names){
		            mtext(trans[i], side=2, adj=0.5, padj=1, line=0, outer=FALSE, cex=0.8, las=1)
			    }
			}
			}
		}
		axis(1, at=round(seq(min(sub$start), max(sub$end), length.out=10)), labels=round(seq(min(sub$start), max(sub$end), length.out=10)), pos=0, cex.axis=1.5)
	}
	if(legend){
		mtext(paste(geneID, unique(featureData(ecs)$strand[rt])), side=3, adj=0.25, padj=1.5, line=0, outer=TRUE, cex=1.5)
		posforlegend <- seq(.7, .9, length.out=numcond)
		for(i in 1:length(color)){
			mtext(names(color[i]), side=3, adj=posforlegend[i], padj=1.5, line=0, outer=TRUE, cex=1.5, col=color[i])
		}
	}else{
	mtext(paste(geneID, unique(featureData(ecs)$strand[rt])), side=3, adj=0.5, padj=1.5, line=0, outer=TRUE, cex=1.5)
	}
}
