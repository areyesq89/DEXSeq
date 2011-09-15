\name{estimatelog2FoldChanges}
     \alias{estimatelog2FoldChanges}
     \title{Fold changes (log2) from the coefficients of a glm.}
     \description{
	This function adds to the featureData of an ExonCountSet object the log fold changes (log2) of the expression levels.  It calculates them from the coefficients of the glm.
     }
     \usage{
     estimatelog2FoldChanges(ecs, fitExpToVar="condition", n.cores=1, quiet=FALSE)
     }
     \arguments{
       \item{ecs}{An ExonCountSet object.}
       \item{fitExpToVar}{A variable contained in the design of the ecs, the expression values will be fitted to this variable using the the formula "count~sample + fitExpToVar * exon"}
       \item{n.cores}{Number of cores to be used to estimate the dispersions.  multicore package must be loaded in order to split the job in several cores.}
       \item{quiet}{If TRUE, no progress report is shown.}
     }
     \seealso{
       graphics, \code{segments}
     }
     \examples{
	\dontrun{data("pasillaExons", package="pasilla")}
	\dontrun{pasillaExons <- estimateSizeFactors(pasillaExons)}
	\dontrun{pasillaExons <- estimateDispersions(pasillaExons)}
	\dontrun{plotDEXSeq(pasillaExons, "FBgn0085442")}
	\dontrun{pasillaExons <- estimatelog2FoldChanges(pasillaExons, fitExpToVar="condition", n.cores=1)}
     }
\keyword{ExonCountSet}
