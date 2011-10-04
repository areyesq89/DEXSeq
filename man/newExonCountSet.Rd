\name{newExonCountSet}
     \alias{newExonCountSet}
     \title{Creates an ExonCountSet object}
     \description{
       This function creates an ExonCountSet object from a matrix or data.frame of read counts.
     }
     \usage{
newExonCountSet(countData, design, geneIDs, exonIDs, exonIntervals=NULL, 
                           transcripts=NULL)
     }
     \arguments{
       \item{countData}{
           A matrix or data frame of count data of non-negative integer values. 
           The rows correspond to counts for each exon, the columns correspond to samples 
           (or experiments). Note that biological replicates should each get their own column, 
           while the counts of technical replicates (i.e., several sequencing ruins/lanes from 
           the same sample) have to be summed up into a single column.}
       \item{design}{
           A factor or data frame with the design annotation (e.g. treatments, or tissue 
           types, or phenotypes, or the like).  The length of the factor (or rows in the data frame) 
           has to be equal to the number of columns of the countData matrix, assigning a 
           condition to each sample. If 'design' is not a factor, it will be converted to one.
        }
	\item{geneIDs}{
           A vector of gene identificators ordered according to its respective row in countData.  
           If the gene "x" has four exons and therefore four columns in countData, then "x" must 
           be four times in the vector. If it is not a factor, it will be converted to one.
        }
	\item{exonIDs}{
           A  character vector of exon identificators ordered according to its respective row in 
           countData.  The identificator names can be repeated between genes but not within genes. 
        }
	\item{exonIntervals}{
           Receives a data frame with exon annotation information. The number of rows in the data 
           needs to be of the same length as the number of rows in countData. The columns names must 
           contain the values "chr", "start", "end", "strand".  It is not extrictly necessary but 
           it must be added in order to use the \code{plotDEXSeq} function.
        }
	\item{transcripts}{
           A character vector of the same length as the rows of the count data containing, for each 
           row in countData, a concatenation of transcriptIDs separated by the character ";".  This 
           means that if an exon is contained in the transcripts "A", "B" and "C", the field of the 
           row corresponding to that exon should contain "A;B;C".
        }
     }
     \value{
	An object of class ExonCountSet.
     }
     \seealso{
       \code{\link{read.HTSeqCounts}}
     }
     \examples{
data("pasillaExons", package="pasilla")
ecs <- newExonCountSet(
   countData=counts(pasillaExons), 
   design=design(pasillaExons), 
   geneIDs=geneIDs(pasillaExons), 
   exonIDs=exonIDs(pasillaExons))
}
\keyword{ExonCountSet}
