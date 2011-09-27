\name{estimatelog2FoldChanges}
   \alias{estimatelog2FoldChanges}
   \title{Fold changes (log2) from the fitted expression values in the glm.}
   \description{
      This function calculates the fold changes (log2) between the different
      conditions. It calculates them from the coefficients of a glm 
      that fits the read counts to a variable of the experimental design specified 
      by the user.
   }
   \usage{
      estimatelog2FoldChanges(ecs, fitExpToVar="condition", nCores=1, quiet=FALSE, file="")
   }
   \arguments{
      \item{ecs}{
         An ExonCountSet object.
      }
      \item{fitExpToVar}{
         A variable contained in the design of the ecs, the expression 
         values will be fitted to this variable using the the formula 
         "count~sample + fitExpToVar * exon".
      }
      \item{nCores}{
         Number of cores to be used to estimate the dispersions.  
         multicore package need to be loaded before in order to split the 
         job in several cores.
      }
      \item{quiet}{
         If TRUE, no progress report is shown.  In case the session is 
         not an interactive and progress report is wanted.
         Change the name of the file.
      }
      \item{file}{
         A file name to write the progress reports. If file="", output
         will be written in the standart output connection.
      }
  }
  \examples{
  \dontrun{data("pasillaExons", package="pasilla")}
  \dontrun{pasillaExons <- estimateSizeFactors(pasillaExons)}
  \dontrun{pasillaExons <- estimateDispersions(pasillaExons)}
  \dontrun{plotDEXSeq(pasillaExons, "FBgn0085442")}
  \dontrun{pasillaExons <- estimatelog2FoldChanges(pasillaExons, fitExpToVar="condition", nCores=1)}
  }
\keyword{ExonCountSet}
