\name{DEXSeqHTML}
     \alias{DEXSeqHTML}
     \title{HTML report writer}
     \description{
       This function makes a report from the results from \code{testForDEU} saved in an ExonCountSet object. It allows an easy way of exploring the results of the test.
     }
     \usage{
     DEXSeqHTML(ecs, geneIDs=NULL, path="DEXSeqReport", file="testForDEU.html", fitExpToVar="condition", FDR=0.1, color=NULL, color.samples=NULL)
     }
     \arguments{
       \item{ecs}{An ExonCountSet object}
       \item{geneIDs}{A character vector of gene identificators to be included in the report. The default is NULL and the default genes included in the report are the ones with significant changes in exon usage below the FDR.}
       \item{path}{A path where to write the HTML file and plots for the visualization of the results of \code{testForDEU}.}
       \item{file}{A file to write the html report}
       \item{fitExpToVar}{A variable contained in the design of the ecs, the expression values will be fitted to this variable}
       \item{FDR}{A false discovery rate for the results.}
       \item{color}{A vector of colors for each of the conditions indicated by fitExpToVar.}
       \item{color.samples}{A vector of colors for each of the samples. If NULL, the colors of each sample will be asigned according to its corresponding condition. Useful to visualize complex experimental designs.}
     }
     \value{
	This function will write an report in the specified path. In this path it will create an html file with the initial report page and a path called "files" to which the plots and other htmls are going to be created.  To see an example please visit \url{http://www.embl.de/~reyes/DEXSeqReport/testForDEU.html}.
     }
     \seealso{
       \code{hwrite}
     }
     \examples{
	\dontrun{data("pasillaExons", package="pasilla")}
	\dontrun{pasillaExons <- estimateSizeFactors( pasillaExons )}
	\dontrun{pasillaExons <- estimateDispersions( pasillaExons )}
	\dontrun{pasillaExons <- testForDEU( pasillaExons )}
	\dontrun{DEXSeqHTML( pasillaExons )}
     }
\keyword{ExonCountSet}
