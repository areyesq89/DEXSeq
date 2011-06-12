\name{DEUresultTable}
\alias{DEUresultTable}
\title{
      Makes a result table from \code{testForDEU}.
}
\description{
      This function makes a table describing the dispersion parameter used in each of the test and the results from \code{testForDEU}.
}
\usage{
DEUresultTable(ecs, foldChange=FALSE)
}
\arguments{
  \item{ecs}{
     An ExonCountSet object.
	}
  \item{foldChange}{
     A logic parameter indicating whether or not the log2 fold change should be included in the report.  It is only available when comparing between two conditions.
	}
}
\value{
	A data frame that summarizes the results from \code{testForDEU}.
}
