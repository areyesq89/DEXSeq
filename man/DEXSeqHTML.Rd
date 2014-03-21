\name{DEXSeqHTML}
\alias{DEXSeqHTML}
\title{DEXSeq HTML report writer}
\description{
  This function generates an HTML report from the results
  of a \code{DEXSeqResults} object. Gives an simple report
  to explore differential exon usage results.
}
\usage{
  DEXSeqHTML(object, genes=NULL, path="DEXSeqReport", file="testForDEU.html", 
    fitExpToVar="condition", FDR=0.1, color=NULL, color.samples=NULL,
    mart="", filter="", attributes="", extraCols=NULL, BPPARAM=MulticoreParam(workers=1))
}
\arguments{
  \item{object}{
    An DEXSeqResults object
  }
  \item{genes}{
    A character vector of gene identifiers 
    to be included in the report. If NULL, 
    the genes included in the report will 
    be the significant hits at the given 
    false discovery rate. See "FDR" below.
  }
  \item{path}{
    A path in the system where to write the 
    report.
  }
  \item{file}{
    The name of the html file.
  }
  \item{fitExpToVar}{
    A variable contained in the design of the ecs; the counts 
    will be fitted to this variable to get the plotting values.
  }
  \item{FDR}{
    A false discovery rate
  }
  \item{color}{
    A vector of colors, one for each of the levels of the 
    values of "fitExpToVar".
  }
  \item{color.samples}{
    A vector of colors for each of the samples. If NULL, the colors 
    of each sample will be asigned according to its corresponding 
    condition. Useful to visualize complex experimental designs.
   }
   \item{mart}{
    object of class Mart, created with the useMart function, 
    with dataset specified
   }
   \item{filter}{ 
     Filters (ONLY ONE) that should be used in the query.  A
     possible list of filters can be retrieved using the function
     listFilters. Please note that the value of this filter will
     always be the geneIDs in the DEXSeqResults object.
   }
   \item{attributes}{
     Attributes you want to retrieve.  A possible list of
     attributes can be retrieved using the biomaRt function
     listAttributes.
   }
   \item{extraCols}{
     A data frame with one or more columns to add to the report.
     For example, additional information about the genes. The data frame
     should be indexed by the gene names of the ExonCountSet object,
     e.g. the rownames of the data frame should correspond to the gene 
     names.
   }
  \item{BPPARAM}{ 
     A 'BiocParallelParam' instance. 
     See \code{?bplapply} for details.
  }
}
\value{
  This function will write an HTML report in the directory specified by 'path'. 
  There, it will create an html file with the initial report page and a 
  directory called "files" in which SVG files with the plots and other html files are
  placed. Different plots with different labels are generated for each gene:
   - counts: the raw data, for each sample
   - fitted expression: the fitted coefficients per compared condition (e.g.: treated, untreated)
   - fitted splicing: as 'expression', but after removing overall gene-level differential expression: this is the view most relevant for the interpretation of DEXSeq results, which are about changes in relative exon usage (i.e.: relative to overall gene expression) 
}
\examples{
	\dontrun{

        }
}
