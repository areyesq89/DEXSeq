\name{DEXSeqHTML}
     \alias{DEXSeqHTML}
     \title{HTML report writer}
     \description{
      This function generates an HTML report from the results from 
      \code{\link{testForDEU}} saved in an ExonCountSet object.
      It uses the information from the function 
      \code{\link{DEUresultTable}} and plotting from \code{\link{plotDEXSeq}}.
      It allows an easy way of exploring the results of the tests.
     }
     \usage{
     DEXSeqHTML(ecs, geneIDs=NULL, path="DEXSeqReport", file="testForDEU.html", fitExpToVar="condition", FDR=0.1, color=NULL, color.samples=NULL)
     }
     \arguments{
       \item{ecs}{An ExonCountSet object}
       \item{geneIDs}{
           A character vector of gene identificators to be included in the report. 
           If leaved as NULL, the genes included in the report will be the significant
           hits at the given false discovery rate. See below "FDR".}
       \item{path}{
           A path in the system where to write the report.
       }
       \item{file}{
           The name of the html file.
       }
       \item{fitExpToVar}{
           A variable contained in the design of the ecs, the counts 
           will be fitted to this variable for the visualization.
       }
       \item{FDR}{A false discovery rate for the result.}
       \item{color}{A vector of colors for each of the levels of the value from fitExpToVar.}
       \item{color.samples}{
           A vector of colors for each of the samples. If NULL, the colors 
           of each sample will be asigned according to its corresponding 
           condition. Useful to visualize complex experimental designs.
       }
     }
     \value{
	   This function will write an HTML report in the specified path. In this 
           path it will create an html file with the initial report page and a 
           path called "files" to which the plots and other htmls are going to 
           be created.  To see an example please visit 
           \url{http://www.embl.de/~reyes/DEXSeqReport/testForDEU.html}.
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
