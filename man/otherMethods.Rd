\name{DEXSeqDataSet-methods}
\docType{methods}
\alias{[,DEXSeqDataSet-method}
\alias{[,DEXSeqDataSet,ANY-method}
\alias{[,DEXSeqDataSet,ANY,ANY-method}
\title{Subset DEXSeqDataSet object}
\description{
  Subsets a DEXSeqDataSet object.
}
\usage{
    \S4method{[}{DEXSeqDataSet}(x, i, j, ..., drop=TRUE)
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
    \item{...}{
        Other arguments passed to lower functions
    }
    \item{drop}{
        Ignored for subsetting DEXSeqDataSet objects.
    }
}