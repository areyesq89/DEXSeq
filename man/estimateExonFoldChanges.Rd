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
     BPPARAM=MulticoreParam(workers=1), maxRowsMF=3000)

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
  \item{maxRowsMF}{
     For the fold change estimation, the size of the model matrix 
     for the fitted glm increases with the number of samples and 
     the number of exons for a specific gene (see the DEXSeq paper
     for details). Since the glm fit for big models can be slow, 
     the maxRowsMF parameter allows to set a threshold on the number
     of rows from the model frame (this will be number of samples x number
     of exons from a gene).  For the genes passing this threshold 
     the exon fold changes won't be calculated.
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
  data(pasillaDEXSeqDataSet, package="pasilla")
  dxd <- estimateSizeFactors( dxd )
  dxd <- estimateDispersions( dxd )
  dxd <- testForDEU( dxd )
  dxd <- estimateExonFoldChanges( dxd )
}

