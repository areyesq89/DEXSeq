\name{estimateExonFoldChanges}
\alias{estimateExonFoldChanges}
\title{Estimates exon usage coefficients from the fitted terms of the GLM.}
\description{
  This function calculates the exon usage coefficients and fold changes 
  (on log2 scale) between the different conditions.
}
\usage{

  estimateExonFoldChanges( object, 
     fitExpToVar = "condition", denominator = "",
     BPPARAM=MulticoreParam(workers=1))

}
\arguments{
  \item{object}{
    An DEXSeqDataSet object.
  }
  \item{fitExpToVar}{
    A variable contained in \code{design(ecs)}. The expression 
    values will be fitted to this variable using the the formula 
    " ~ sample + fitExpToVar * exon".
  }
  \item{denominator}{
    A value of the sample annotation (e.g. condition) to use as a denominator 
    in the log2 fold change. As a default, the function will take 
    the annotation of the first sample.
  }
  \item{BPPARAM}{
     A "BiocParallelParam" instance. 
     See \code{?bplapply} for details.
  }
}
\details{
  Exon usage coefficients are calculated by fitting
  a GLM from the joint data of all exons of the same gene.
  The model frame can be found in the slot object@modelFrameBM
  of a \code{DEXSeqDataSet} object.  The model 
  '~ sample + fitExpToVar * exon' is fitted.  
  The resulted coefficients are arranged and reformatted 
  in order to remove gene expression effects (absorbed by the
  'sample' variable in the formula), leaving only exon 
  usage effects for each individual exon in each level of
  'fitExpToVar'. These values are used by the function \code{plotDEXSeq}. 
}
\examples{

}

