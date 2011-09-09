\name{plotDEXSeq}
     \alias{plotDEXSeq}
     \title{Visualization of the coefficients of the glm or the normalized counts.}
     \description{
	The function provides a tool to visualize the gene annotation information, read count data and the results of the test in \code{testForDEU}. It is made on graphics.
     }
     \usage{
     plotDEXSeq(ecs, geneID, FDR=0.1, fitExpToVar="condition", norCounts=FALSE, expression=TRUE, splicing=FALSE, displayTranscripts=FALSE, names=FALSE, legend=FALSE, color=NULL, color.samples=NULL, ...)
     }
     \arguments{
       \item{ecs}{An ExonCountSet object.}
       \item{geneID}{geneID to visualize.}
       \item{FDR}{A false discovery rate threshold used to indicate the significant exons.}
       \item{fitExpToVar}{A variable contained in the design of the ecs, the expression values will be fitted to this variable}
       \item{norCounts}{If TRUE, provides a plot of the normalized counts by the size factors.}
       \item{expression}{If TRUE, the function plots of the fitted EXPRESSION estimates from the glm regression.}
       \item{splicing}{If TRUE, the function plots of the fitted SPLICING estimates from the glm regression.}
       \item{displayTranscripts}{Logic argument which indicates if the annotated transcripts want to be plotted or not, the upper limit is 42 transcripts due to graphics limits.}
       \item{names}{Logical indicating if the names of the transcriptIDs should to be written.}
       \item{legend}{Logical indicating if the legend of colors and conditions should be added.}
       \item{color}{A vector of colors for each of the conditions indicated by fitExpToVar.}
       \item{color.samples}{A vector of colors for each of the samples. If NULL, the colors of each sample will be asigned according to its corresponding condition. Useful to visualize complex experimental designs.}
       \item{...}{Further graphical parameters (from par)}
     }
     \seealso{
       graphics, \code{segments}
     }
     \examples{
	\dontrun{data("pasillaExons", package="pasilla")}
	\dontrun{pasillaExons <- estimateSizeFactors(pasillaExons)}
	\dontrun{pasillaExons <- estimateDispersions(pasillaExons)}
	\dontrun{plotDEXSeq(pasillaExons, "FBgn0085442")}
     }
\keyword{ExonCountSet}
