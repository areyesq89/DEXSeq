\name{estimateDispersions}
\docType{methods}
\alias{estimateDispersions,ExonCountSet-method}
\alias{estimateDispersions}
\title{Estimate exon dispersions}
\description{
      This function estimates the Cox-Reid profile log likelihood 
      dispersion estimates for each exon of each of the genes of 
      the ExonCountSet object. It stores this values in 
      \code{fData(ecs)$dispersionBeforeSharing}.  
}
\usage{
\S4method{estimateDispersions}{ExonCountSet}( cds, formula=count ~ sample + condition*exon, initialGuess=.01, nCores=1, minCount=10, maxExon=70, quiet=FALSE, file="")
}
\arguments{
  \item{cds}{
     An ExonCountSet object.
   }
  \item{formula}{
     Formula used in the glm to calculate the Cox-Reid dispersion values.  
     The factors on the formula must be present in the design 
     columns of the ExonCountSet object.
   }
  \item{initialGuess}{
     An arbitrary initial guess for the dispersion values to initiate the optimization.
   }
  \item{nCores}{
     Number of cores to be used to estimate the dispersions.  
     Multicore package must be loaded in order to split the 
     job in several cores.
   }
  \item{minCount}{
     Minimum number of counts on an exon for it to be considered in the tests. 
     This significantly increases the speed of the dispersion 
     estimations and testing for differential exon usage. 
     This is supported by the fact that small count exons are less likely
     of being called significant, so it should not affect the results.
   }
  \item{maxExon}{
     Genes with more exons than this value will be discarded from 
     the analysis. This is a speed issue. Currently, time of dispersion 
     estimations and testing for differential exon usage
     increases with number of exons.
   }
   \item{quiet}{
     If TRUE, no progress report is shown.  In case the session is 
     not an interactive session and progress report is wanted, include
     a file name in the parameter "file".
   }
   \item{file}{
     A file name to write the progress reports. If file="", output
     will be written in the standart output connection.
   }
}
\value{
     An ExonCountSet object with dispersion \code{featureData(cds)$dispersion_CR_est}) 
     parameters filled).
}
\examples{
     \dontrun{data("pasillaExons", package="pasilla")}
     \dontrun{pasillaExons <- estimateSizeFactors( pasillaExons )}
     \dontrun{pasillaExons <- estimateDispersions( pasillaExons )}
}
