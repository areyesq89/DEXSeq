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
     BPPARAM=SerialParam(), maxRowsMF=2400,
	 independentFiltering=TRUE, filter)

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
     For the exon fold change estimation, the size of the model matrix 
     for the fitted glm increases with the number of samples and 
     the number of exons for a specific gene (see the DEXSeq paper
     for details). Since the glm fit for big models is very slow, 
     the maxRowsMF parameter allows to set a threshold on the number
     of rows from the model matrix (the number of rows of the model
     matrix will be number of samples times number of exons of a
     gene). For all genes passing this threshold, the exon fold
     changes will be estimated by fitting a slightly different
     but equivalent model that considers, for each exon, the
     counts from the specific exon and the sum of the rest of exons
     of the same gene.
  }
  \item{independentFiltering}{
     Logical indicating whether independent filtering should be 
	 applied automatically. 
  }
  \item{filter}{
     A vector of filter statistics over which the independent
	 filtering will be optimized. The default is the normalized
	 exon means.
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

