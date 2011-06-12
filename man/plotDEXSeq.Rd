\name{plotDEXSeq}
     \alias{plotDEXSeq}
     \title{Visualization of the coefficients of the glm and/or the normalized counts.}
     \description{
	The function provides a tool to visualize the gene annotation information, read count data and the results of the test in \code{testForDEU}. It is made on graphics.
     }
     \usage{
     plotDEXSeq(ecs, geneID, FDR=NULL, coefficients=TRUE, padj=NULL, norCounts=FALSE, displayTranscripts=FALSE, names=FALSE, color=NULL, ...)
     }
     \arguments{
       \item{ecs}{An ExonCountSet object.}
       \item{geneID}{geneID to visualize.}
       \item{FDR}{A false discovery rate threshold used to indicate the significant exons.}
       \item{coefficients}{If TRUE, provides a plot of the log2 coefficient estimates from the glm regression, which are an estimate of the mean of the number of reads falling in the exon.}
       \item{norCounts}{If TRUE, provides a plot of the normalized counts by the SizeFactor.}
       \item{expression}{Logical argument, if TRUE, the effect of the condition in changes in the expression in the gene is included in the plot. IF FALSE, pure splicing effect is plotted.}
       \item{displayTranscripts}{Logic argument which indicates if the annotated transcripts want to be plotted or not, the upper limit is 45 transcripts.}
       \item{color}{A vector of colors indicating a color foreach condition.}
       \item{padj}{A vector of adjusted pvalues, if NULL they will be calculated using the function \code{testGeneForDEU}}
       \item{legend}{Logical. A legend indicating the colors for each condition will be added.}
       \item{names}{Logical whether or not indicating if the names of the transcript identificators want to be written.}
       \item{...}{Further graphical parameters (from ‘par’)}
     }
     \seealso{
       graphics, \code{segments}
     }
     \examples{
	makeExampleExonCountSet(ecs)
	ecs <- estimateSizeFactors(ecs)
	ecs <- estimateDispersions(ecs)
	plotDEXSeq(ecs, "SG10020")
     }
\keyword{ExonCountSet}
