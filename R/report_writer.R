DEXSeqHTML <- function(ecs, geneIDs=NULL, path="DEXSeqReport", file="testForDEU.html", FDR=0.05, color=NULL){
        stopifnot(is(ecs, "ExonCountSet"))
#       assertOneWay(ecs)
	results<-DEUresultTable(ecs)
	results$dispersion <- round(results$dispersion, 4)
	results$pvalue <- round(results$pvalue, 3)
	results$padjust <- round(results$padjust, 3)
	if(!is.null(results$log2change)){
		results$log2change <- round(results$log2change, 3)
	}
	
	results[which(is.na(results$pvalue)),][,c("pvalue","padjust")]=1
	rownames(results) <- NULL
	numcond <- length(unique(design(ecs)))
	sortabletag <- hmakeTag(tag="script", src=paste(system.file(package="DEXSeq"), "/sorttable.js", sep=''))
	if(is.null(color)){
		color<-rgb(colorRamp(c("#D7191C", "#FFFFBF", "#2B83BA"))(seq(0, 1, length.out=numcond)), max=255, alpha=175)
	}
	dir.create(path)
	dir.create(paste(path, "/", "files", sep=""))
	if(is.null(geneIDs)){
		gns <- as.character(unique(results$geneID[which(results$padjust <= FDR)]))
	}else{
		gns <- geneIDs
	}
	if(length(gns)==0){ stop("There are no significant results in the test... nothing to report")}
	p<-openPage(paste(path, "/", file, sep=""))
	hwrite(sortabletag, p)
	hwrite('DEXSeq differential exon usage test', p, heading=1)
	hwrite('Experimental design', p, heading=2)
	cond<-as.matrix(design(ecs))
	colnames(cond)<-ecs@designColumns
	hwrite(cond, p)
	hwrite('testForDEU result table', p, heading=2)
	ptowrite <- paste(path, "/", "files", "/", sep="")

	######### prepare colors for table results
	
	m2col <- colorRampPalette(c("#FF706B", "#FEE08B", "white"), space="rgb")(5)
	matcol <- matrix(rep(m2col[5], ncol(results)*nrow(results)), nrow(results))
	matcol[which(results$padjust <= 0.25),] <- m2col[4]
	matcol[which(results$padjust <= 0.1),] <- m2col[3]
	matcol[which(results$padjust <= 0.05),] <- m2col[2]
	matcol[which(results$padjust <= 0.01),] <- m2col[1]
	legend <- hwrite(c("<= 0.01", "<= 0.05", "<= 0.1", "<= 0.25", "> 0.25"), bgcolor=m2col)



	for( gene in gns){
		back <- hwrite("back", link=paste("../", file, sep=""))
		otherlinks <- hwrite(c("counts", "expression", "splicing", "transcripts", "results"), link=c(paste(gene, "counts.html", sep=""), paste(gene, "expression.html", sep=""), paste(gene, "splicing.html", sep=""), paste(gene, "transcripts.html", sep=""), paste(gene, "results.html", sep="")), table=FALSE)
		
		subres <- results[as.character(results$geneID) %in% as.character(gene),]
		submatcol <- matcol[as.character(results$geneID) %in% as.character(gene),]
		rownames(subres) <- NULL
		genpage <- openPage(paste(ptowrite, gene, "results.html", sep=""))
		hwrite(sortabletag, genpage)
		hwrite(c(back, otherlinks), table=TRUE, border=0, genpage)
		hwrite(legend, p=genpage)
		hwrite(as.matrix(subres), bgcolor=submatcol, table.class="sortable", style='margin:16px; border:0px solid black; border-width:0px; width:200px', table=TRUE, page=genpage)
		close(genpage, splash=TRUE)

		

		genpage <- openPage(paste(ptowrite, gene, "expression", ".html", sep=""))
		hwrite(c(back, otherlinks), table=TRUE, border=0, genpage)
		svg(paste(ptowrite, gene, "expressionsvg.svg", sep=""), height=7, width=12, pointsize=14)
		plotDEXSeq(ecs, gene, FDR, lwd=2, legend=TRUE, color=color, cex.axis=1.5)
		dev.off()
		hwrite(hmakeTag("iframe", src=paste(gene, "expressionsvg.svg", sep=""), width=1200, height=700, border=0), p=genpage)
		close(genpage, splash=TRUE)

		genpage <- openPage(paste(ptowrite, gene, "counts.html", sep=""))
		hwrite(c(back, otherlinks), table=TRUE, border=0, genpage)
		svg(paste(ptowrite, "/", gene, "countssvg.svg", sep=""), height=7, width=12, pointsize=14)
		plotDEXSeq(ecs, gene, FDR, coefficients=FALSE, norCounts=TRUE, lwd=2, legend=TRUE, color=color, cex.axis=1.5)
		dev.off()
		hwrite(hmakeTag("iframe", src=paste(gene, "countssvg.svg", sep=""), width=1200, height=700, border=0), genpage)
		close(genpage, splash=TRUE)


		genpage <- openPage(paste(ptowrite, gene, "splicing.html", sep=""))
		hwrite(c(back, otherlinks), table=TRUE, border=0, genpage)
		svg(paste(ptowrite, gene, "splicingsvg.svg", sep=""), height=7, width=12, pointsize=14)
		plotDEXSeq(ecs, gene, FDR, expression=FALSE, lwd=2, legend=TRUE, color=color, cex.axis=1.5)
		dev.off()
		hwrite(hmakeTag("iframe", src=paste(gene, "splicingsvg.svg", sep=""), width=1200, height=700, border=0), genpage)
		close(genpage, splash=TRUE)
		
		if(!is.null(featureData(ecs)$transcripts)){
			transcripts <- sapply(featureData(ecs)$transcripts[featureData(ecs)$geneID %in% gene], function(x){strsplit(x, ";")})
			trans <- Reduce(union, transcripts)
			if(length(trans) > 10){    ############ if there are more than 10 transcripts, increase the size of the plotting region
				h <- (7+(length(trans)*0.3))
			}else{
				h <- 7
			}
			if(sum(as.character(results$geneID) %in% as.character(gene)) > 35){   ############# if there are more than 35 exons, increase the size of the plotting region
				h <- h + (sum(as.character(results$geneID) %in% as.character(gene))*0.05)
			}
			genpage <- openPage(paste(ptowrite, gene, "transcripts", ".html", sep=""))
			hwrite(c(back, otherlinks), table=TRUE, border=0, genpage)
			svg(paste(ptowrite, "/", gene, "transcriptssvg.svg", sep=""), height=h, width=12, pointsize=14)
			plotDEXSeq(ecs, gene, FDR, displayTranscripts=TRUE, lwd=2,  legend=TRUE, color=color, cex.axis=1.5)
			dev.off()
			hwrite(hmakeTag("iframe", src=paste(gene, "transcriptssvg.svg", sep=""), width=1200, height=h*100, border=0), genpage)
			close(genpage, splash=TRUE)
		}
	}
	
	
	results <- results[as.character(results$geneID) %in% gns,]
	genetable <- cbind( 
				 geneID=unique(as.character(results$geneID)),
				 do.call(rbind, 
						lapply(unique(as.character(results$geneID)), function(gene){
							vec <- as.character(geneIDs(ecs)) %in% gene 
							data.frame(chr=unique(featureData(ecs)$chr[vec]), start=min(featureData(ecs)$start[vec]), end=max(featureData(ecs)$end[vec]))
							}
						)
				),
				total_exons = rle(as.character(results$geneID))$lengths,
				exon_changes = sapply(unique(as.character(results$geneID)), function(gene){vec <- as.character(results$geneID) %in% gene; sum(results$padjust[vec] <= FDR)})
			  )
	
	genetable$geneID <- sapply(as.character(genetable$geneID), function(m){hwrite(m, link=paste("files/", m, "expression.html", sep=""))})
	rownames(genetable) <- NULL
	hwrite(genetable, page=p, table=TRUE, table.class="sortable", style='margin:16px; border:0px solid black; border-width:0px; width:200px') 
	close(p, splash=TRUE)
}
