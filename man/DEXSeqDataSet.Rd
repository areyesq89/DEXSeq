\docType{class}
\name{DEXSeqDataSet-class}
\alias{DEXSeqDataSet}
\alias{DEXSeqDataSet-class}
\alias{DEXSeqDataSetFromHTSeq}
\alias{DEXSeqDataSetFromSE}
\title{DEXSeqDataSet object and constructors}
\usage{

  DEXSeqDataSet( countData, sampleData, 
    design= ~ sample + exon + condition:exon, 
    featureID, groupID, featureRanges=NULL, 
    transcripts=NULL)

  DEXSeqDataSetFromHTSeq( 
    countfiles, sampleData, 
    design= ~ sample + exon + condition:exon, 
    flattenedfile=NULL )

  DEXSeqDataSetFromSE( SE, 
    design= ~ sample + exon + condition:exon )

}
\arguments{
  \item{countData}{
    A matrix of count data of non-negative integer values. 
    The rows correspond to counts for each exon counting bin, 
    the columns correspond to samples. Note that biological 
    replicates should each get their own column, while the 
    counts of technical replicates (i.e., several sequencing 
    runs/lanes from the same sample) should be summed up 
    into a single column
  }
  \item{countfiles}{
    A character vector containing the path to the files that
    were originated with the script 'dexseq_count.py'.
  }
  \item{sampleData}{
    A \code{data.frame} with the annotation (e.g. treatments, or tissue 
    types, or phenotypes, or the like).  The number of rows in 
    the data frame must to be equal to the number of columns of 
    the countData matrix, assigning the annotation of each sample.
  }
  \item{design}{
    A formula which specifies the design of the experiment. It must
    specify an interaction term between a variable from the sampleData
    columns with the 'exon' variable. By default, the design 
    will be '~ sample + exon + condition:exon'. This formula indicates
    the contrast between 'condition' and exon', i.e. differences in exon 
    usage due to changes in the 'condition' variable.  See the vignette 
    for more examples of other designs.
  }
  \item{featureID}{
    A character vector of counting regions identifiers ordered according 
    to the rows in countData.  The identifiers names can be repeated
    between groups but not within groups. 
  }
  \item{groupID}{
    A vector of group identifiers ordered according to its
    respective row in countData. It must reflect the sets
    of counting regions belonging to the same group, for 
    example, exon bins in belonging to the same gene should
    have the same group identifier. 
  }
  \item{featureRanges}{
    Optional. \code{GRanges} or \code{GRangesList} describing
    the genomic coordinates of each of the rows of countData.
  }
  \item{transcripts}{
    Optional. A \code{list} of the same length as the number of
    columns of countData. Each element of the list should contain
    the transcript identifiers of the 
  }
  \item{flattenedfile}{
    A character vector containing the path to the flattened annotation 
    file that was originated with the script 'dexseq_prepare_annotation.py'.
  }
  \item{SE}{
    A \code{SummarizedExperiments} object, originated using the function
    \code{SummarizeOverlaps}. See examples for more details.
  }  
}

\value{
  A DEXSeqDataSet object.
}

\seealso{
  \code{DEXSeqDataSetFromHTSeq}
  \code{DEXSeqDataSetFromSE}
}

\description{
  The \code{DEXSeqDataSet} is a subclass of \code{DESeqDataSet},
  specifically designed to adapt the \code{DESeqDataSet} to test
  for differences in exon usage. 
}
\examples{

\dontrun{
  #######################################
  ### From the output of the          ### 
  ### acconpaning python scripts      ### 
  #######################################

  inDir = system.file("extdata", package="pasilla", mustWork=TRUE)
  flattenedfile = file.path(inDir, "Dmel.BDGP5.25.62.DEXSeq.chr.gff")
  sampleData = data.frame(
      condition = c( rep("treated", 3), rep("untreated", 4) ),
      type = c("single", "paired", "paired", "single", "single", "paired", "paired") )
     
  countFiles <- list.files(inDir, pattern="fb.txt")
  rownames( sampleData ) <- countFiles

  DEXSeqDataSetFromHTSeq(
    countfiles=file.path( inDir, countFiles ),
    sampleData = sampleData,
    design = ~ sample + exon + type:exon + condition:exon,
    flattenedfile=flattenedfile )


  #######################################
  ### From GRanges derived objects    ###
  #######################################

  library(GenomicRanges)
  library(GenomicFeatures)
  library(GenomicAlignments)

  hse <- makeTranscriptDbFromBiomart( 
    biomart="ensembl", dataset="hsapiens_gene_ensembl" )

  bamDir <- system.file( 
    "extdata", package="parathyroidSE", mustWork=TRUE )
  fls <- list.files( bamDir, pattern="bam$", full=TRUE )
  
  bamlst <- BamFileList( 
    fls, index=character(), 
    yieldSize=100000, obeyQname=TRUE )

  SE <- summarizeOverlaps( exonicParts, bamlst, 
    mode="Union", singleEnd=FALSE,
    ignore.strand=TRUE, inter.feature=FALSE, fragments=TRUE )

  colData(SE)$condition <- c("A", "A", "B")

  DEXSeqDataSetFromSE( SE, 
    design= ~ sample + exon + condition:exon )
  
   

  #######################################
  ### From elementary data structures ###
  #######################################
  countData <- matrix( rpois(10000, 100), nrow=1000 )
  sampleData <- data.frame(
      condition=rep( c("untreated", "treated"), each=5 ) )
  design <- formula( ~ sample + exon + condition:exon )
  groupID <- rep(
      sprintf("gene%d", 1:10),
      each= 100 )
  featureID <- rep(
      sprintf("exon%d", 1:10),
      times= 100 )
  DEXSeqDataSet( countData, sampleData, design,
              featureID, groupID )

}

}
