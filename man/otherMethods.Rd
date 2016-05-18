\name{DEXSeqDataSet-methods}
\docType{methods}
\alias{[,DEXSeqDataSet-method}
\alias{[,DEXSeqDataSet,ANY-method}
\alias{[,DEXSeqDataSet,ANY,ANY-method}
\alias{[,DEXSeqDataSet,ANY,ANY,ANY-method}
\alias{colData<-,DEXSeqDataSet,DataFrame-method}
\alias{$<-,DEXSeqDataSet-method}
\title{Subset and replacements of DEXSeqDataSet slots}
\description{
  Subset and replacements of DEXSeqDataSet slots.
}
\usage{
    \S4method{[}{DEXSeqDataSet}(x, i, j, ..., drop=TRUE)
    \S4method{$}{DEXSeqDataSet}(x, name) <- value
}
\arguments{
    \item{x}{
        A DEXSeqDataSet object.
    }
    \item{i}{
        Indices specifying elements to extract.
    }
    \item{j}{
        Indices specifying elements to extract.
    }
    \item{name}{ A symbol representing the name of a column of
       \code{colData}. }
    \item{value}{An object of a class specified in the S4 method signature.}
    \item{...}{
        Other arguments passed to lower functions
    }
    \item{drop}{
        Ignored for subsetting DEXSeqDataSet objects.
    }
}
