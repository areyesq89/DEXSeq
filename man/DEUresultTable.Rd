\name{DEUresultTable}
\alias{DEUresultTable}
\title{
   Makes a result table from \code{\link{testForDEU}}.
}
\description{
   This function returns a data frame with the summary of the results the analysis workflow.  
   It access the \code{\link{fData}} slots with information of the dispersion estimates obtained
   from the function \code{\link{fitDispersionFunction}}, the p values, and adjusted p values
   obtained from the function \code{\link{testForDEU}}, and log2 fold changes 
   from the function \code{\link{estimatelog2FoldChanges}}.
}
\usage{
   DEUresultTable(ecs)
}
\arguments{
   \item{ecs}{
      An ExonCountSet object.
   }
}
\value{
   A data frame with a summary of the analysis workflow.
}
