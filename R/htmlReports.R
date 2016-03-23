DEXSeqHTML <- function(object, genes=NULL, path="DEXSeqReport", file="testForDEU.html", fitExpToVar="condition", FDR=0.1, color=NULL, color.samples=NULL, mart="", filter="", attributes="", extraCols=NULL, BPPARAM=MulticoreParam(workers=1) )
{
   stopifnot( is( object, "DEXSeqResults" ) )
   if(!fitExpToVar %in% colnames( object@modelFrameBM )){
       stop(sprintf("The value of the parameter fitExpToVar, '%s', is not a variable from the annotation of the samples. Please specify a column name of the colData slot of the DEXSeqDataSet object.", fitExpToVar))
   }
   ######## GET THE RESULT TABLE READY ##########
   genomicData <- as.data.frame( object$genomicData )
   results <- data.frame( object[, c( "groupID", "featureID", "exonBaseMean", "dispersion", "pvalue", "padj" )], stringsAsFactors=TRUE)
   results <- cbind( results, genomicData )

   dexseqR <- mcols( object )$type == "DEXSeq results"

   if(sum(dexseqR, na.rm=TRUE) > 0){
      results <-
          cbind(
              results,
              round(
                  as.data.frame( object[,which(dexseqR)] ), 3
                  )
              )
   }

   rownames(results) <- NULL
   sampleData <- object@sampleData
   
   numcond <- length(unique(sampleData[[fitExpToVar]]))
   if(is.null(color)){
      color<-rgb(colorRamp(c("#D7191C", "#FFFFBF", "#2B83BA"))(seq(0, 1, length.out=numcond)), maxColorValue=255, alpha=175)
   }
   names(color) <- sort(levels(sampleData[[fitExpToVar]]))
   dir.create(file.path(path, "files"), recursive=TRUE)

   if(is.null(genes)){
      gns <- as.character(unique(results$groupID[which(results$padj < FDR)]))
   }else{
      gns <- genes
   }
      
   results[,c("dispersion", "pvalue", "padj")] <- round(results[,c("dispersion", "pvalue", "padj")], 3)
   
   if(!all(gns %in% object$groupID)){
      stop("The geneIDs provided are not in the ecs object")}
   if(length(gns)==0){ 
      stop("There are no significant results in the test... nothing to report")}

   p<-openPage(file.path(path, file))
   hwrite('DEXSeq differential exon usage test', p, heading=1)
   hwrite('Experimental design', p, heading=2)
   cond<-as.matrix( as.data.frame( sampleData[,!colnames(sampleData) %in% "sizeFactor"] ) )
   rownames(cond) <- NULL
   condcolor <- matrix(rep("white", nrow(cond)*ncol(cond)), nrow(cond))
   condcolor[,which(colnames(cond) %in% fitExpToVar)] <- color[as.character(sampleData[[fitExpToVar]])]
   if(!is.null(color.samples)){
      condcolor[,1] <- color.samples}
   hwrite(cond, bgcolor=condcolor, p)

   formulas <- mcols(object)[colnames(object) == "pvalue","description"]
   formulas <- sapply( strsplit(formulas, "vs|p-value:" ), "[", c(2, 3))
   formulas <- as.vector( gsub("'", "", formulas) )
   
   hwrite(paste("\n\nformulaDispersion = ", formulas[1], sep=""), p, heading=3)
   hwrite(paste("\nformula0 = ", formulas[2], sep=""), p, heading=3)
   hwrite(paste("\nformula1 = ", formulas[1], sep=""), p, heading=3)
   hwrite('testForDEU result table', p, heading=2)
   ptowrite <- paste0(path, "/files/")
   ######### prepare colors for table results
	
   m2col <- colorRampPalette(c("#FF706B", "#FEE08B", "white"), space="rgb")(5)
   matcol <- matrix(rep(m2col[5], ncol(results)*nrow(results)), nrow(results))
   ######### COLOR LEGEND = I BET THERE IS A SMARTER WAY OF DOING THIS:
   j <- 4
   for(i in c(0.25, 0.1, 0.05, 0.01)){
      matcol[which(results$padjust <= i),] <- m2col[j]
      j <- j-1
   }
   legend <- hwrite(c("<= 0.01", "<= 0.05", "<= 0.1", "<= 0.25", "> 0.25"), bgcolor=m2col)

  makePagesForGene <- function(gene){
#      print( gene )
      back <- hwrite("back", link=file.path("..", file))
      nameforlinks <- sapply(strsplit(gene, "\\+"), "[[", 1)
      otherlinks <- hwrite(c("counts", "expression", "splicing", "transcripts", "results"), link=c(paste(nameforlinks, "counts.html", sep=""), paste(nameforlinks, "expression.html", sep=""), paste(nameforlinks, "splicing.html", sep=""), paste(nameforlinks, "transcripts.html", sep=""), paste(nameforlinks, "results.html", sep="")), table=FALSE)
      loc <- as.character(results$groupID) %in% as.character(gene)
      ### this makes the page where to explore the pvalues ###
      subres <- results[loc,]
      submatcol <- matcol[loc,]
      rownames(subres) <- NULL
      genpage <- openPage(paste(ptowrite, nameforlinks, "results.html", sep=""))
      hwrite(c(back, otherlinks), table=TRUE, border=0, genpage)
      hwrite(legend, page=genpage)
      hwrite(as.matrix(subres), bgcolor=submatcol, table.class="sortable", style='margin:16px; border:0px solid black; border-width:0px; width:200px', table=TRUE, page=genpage)
      close(genpage, splash=TRUE)
      ### MAKE THE PLOT HTML PAGES FOR expression, counts and splicing
      makePlotPage( object, ptowrite=ptowrite, gene=gene, whichtag="expression", links=c(back, otherlinks), color=color, color.samples=color.samples, FDR=FDR, fitExpToVar=fitExpToVar, width=1200, height=700, h=7)
      makePlotPage( object, ptowrite=ptowrite, gene=gene, whichtag="counts", links=c(back, otherlinks), color=color, color.samples=color.samples, FDR=FDR, fitExpToVar=fitExpToVar, width=1200, height=700, h=7)
      makePlotPage( object, ptowrite=ptowrite, gene=gene, whichtag="splicing", links=c(back, otherlinks), color=color, color.samples=color.samples, FDR=FDR, fitExpToVar=fitExpToVar, width=1200, height=700, h=7)

      transcripts <- object$transcripts[loc]
      if(length( unlist( transcripts ) ) > 0 ) {
         trans <- Reduce(union, transcripts)
         h <- ifelse(length(trans) > 10, 7+(length(trans)*0.3), 7)
      if(sum(loc) > 30){   ############# if there are more than 30 exons, increase the size of the plotting region
         h <- h + (sum(loc)*0.1)
      }
      makePlotPage( object, ptowrite=ptowrite, gene=gene, whichtag=c("expression", "transcripts"), links=c(back, otherlinks), color=color, color.samples=color.samples, FDR=FDR, fitExpToVar=fitExpToVar, width=1200, height=h*100, h=h)
      }
      return()
   }


   bplapply( gns, makePagesForGene, BPPARAM=BPPARAM )
   
   results <- results[as.character(results$groupID) %in% gns,]

   splitCols <- split( seq_len(nrow( results ) ), results$groupID )
   
   genetable <- lapply( splitCols, function(x){
       data.frame(
           chr=unique( results$seqnames[x] ),
           start=min( results$start[x] ),
           end=max( results$end[x] ),
           total_exons = length(x),
           exon_changes = sum( results$padj[x] < FDR, na.rm=TRUE) )
   })
   genetable <- do.call(rbind, genetable)
   genetable <- cbind( geneID=rownames(genetable), genetable )

   if(class(mart) == "Mart"){
      if(attributes(mart)$dataset != ""){
      forvalues <- strsplit(as.character(genetable$geneID), "\\+")
      names(forvalues) <- genetable$geneID
      if(length(filter) > 1){
         warning("length(filter) > 2, only first element will be taken")
         filter <- filter[1]
      }
      extra <- getBM(attributes=c(filter, attributes), filters=filter, values=forvalues, mart=mart)
      fromart <- lapply(genetable$geneID, function(x){
         sep <- do.call(c, strsplit(as.character(x), "\\+"))
         extra[which(extra[,filter] %in% sep),]
      })

      extra <- sapply(attributes, 
         function(r){
           unlist(
              lapply(fromart, 
                 function(x){
                    paste(x[,r], collapse=" ")
                  }
              )
           )
         }
      )
      genetable <- cbind(geneID=genetable$geneID, extra, genetable[,2:length(genetable)])
      }else{
         warning("No dataset in biomart specified")
      }
   }else if( mart != ""){
      warning("Please provide a Mart class object for parameter mart")
   }  

   if( !is.null( extraCols )){
      genetable <- cbind(extraCols[match(genetable$geneID, rownames(extraCols)),], genetable)
   }
   
	
   genetable$geneID <- sapply(as.character(genetable$geneID), function(m){w <- strsplit(m, "\\+");ns <- sapply(w, "[[", 1);hwrite(paste(unlist(w), collapse=" "), link=paste("files/", ns, "expression.html", sep=""))})
   rownames(genetable) <- NULL
   hwrite(genetable, page=p, table=TRUE, table.class="table-layout:fixed", style='margin:16px; border:0px solid black; border-width:1px; width:20%')
   close(p, splash=TRUE)
}


makePlotPage <- function(object, ptowrite, gene, whichtag, links, color, color.samples, FDR, fitExpToVar, width, height, h){
   allopts <- c("expression", "splicing", "counts", "transcripts")
   opts <- allopts %in% whichtag
   onlytag <- allopts[max(which(opts))]
   pagename <- sapply(strsplit(as.character(gene), "\\+"), "[[", 1)
   genpage <- openPage(paste(ptowrite, pagename, onlytag, ".html", sep=""))
   hwrite(links, table=TRUE, border=0, genpage)
   svg(paste(ptowrite, pagename, onlytag, ".svg", sep=""), height=h, width=12, pointsize=14)
   plotDEXSeq(object, geneID=gene, FDR=FDR, lwd=2, expression=opts[1], splicing=opts[2], norCounts=opts[3], displayTranscripts=opts[4], fitExpToVar=fitExpToVar, legend=TRUE, color=color, color.samples=color.samples, cex.axis=1.5)
   dev.off()
   hwrite(hmakeTag("iframe", src=paste(pagename, onlytag, ".svg", sep=""), width=width, height=height, border=0), page=genpage )
   close(genpage, splash=TRUE)
}
     
