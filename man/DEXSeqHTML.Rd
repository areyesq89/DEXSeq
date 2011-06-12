\name{DEXSeqHTML}
     \alias{DEXSeqHTML}
     \title{HTML report writer}
     \description{
       This function makes a report from the results from \code{testForDEU} saved in an ExonCountSet object. It allows an easy way of exploring the results of the test.
     }
     \usage{
     DEXSeqHTML(ecs, geneIDs=NULL, path = "DEXSeqReport", file = "testForDEU.html", FDR = 0.05, foldChange = FALSE, color = NULL)
     }
     \arguments{
       \item{ecs}{An ExonCountSet object}
       \item{geneIDs}{A character vectors of gene identificators to be included in the report. The default is NULL, and the default genes included in the report are the ones with significant changes in exon usage below the FDR.}
       \item{path}{A path where to write the HTML file and plots for the visalization of the results of \code{testForDEU}.}
       \item{file}{A file to write the html report}
       \item{FDR}{A false discovery rate for the results.}
       \item{color}{Vector of colors to be used in the plot for the different conditions.}
     }
     \value{
	This function will write an report in the specified path. In this path it will create an html file with the initial report page and a path called "files" to which the plots and other htmls are going to be created.  To see an example please visit \url{http://www.embl.de/~reyes/DEXSeqReport/testForDEU.html}.
     }
     \seealso{
       \code{hwrite}
     }
     \examples{
	ecs <- makeExampleExonCountSet()
	ecs <- estimateSizeFactors(ecs)
	ecs <- estimateDispersions(ecs)
	ecs <- testForDEU(ecs)
	DEXSeqHTML(ecs)
     }
\keyword{ExonCountSet}
