\name{calcNormFactors}
\docType{methods}
\alias{calcNormFactors,TCC-method}
\alias{calcNormFactors,DGEList-method}
\alias{calcNormFactors}
\title{Calculate normalization factors}
\description{
This function calculates normalization factors using a specified 
multi-step normalization method from a \link{TCC-class} object.
The procedure can generally be described as the 
\eqn{STEP1-(STEP2-STEP3)_{n}}{STEP1-(STEP2-STEP3)n} pipeline.
}
\usage{
\S4method{calcNormFactors}{TCC}(tcc, norm.method = NULL, test.method = NULL,
                iteration = TRUE,  FDR = NULL, floorPDEG = 0.05, 
                increment = FALSE, ...)
}
\arguments{
    \item{tcc}{\link{TCC-class} object.}
    \item{norm.method}{character specifying normalization method that
                       is used in both the STEP1 and STEP3. Possible 
                       values are \code{"tmm"} for the TMM normalization 
                       method implemented in the edgeR package,
                       \code{"edger"} (same as \code{"tmm"}), 
                       and \code{"deseq"} for the method implemented in 
                       the DESeq package. The default is \code{"tmm"}
                       when analyzing the count data with multiple replicates 
                       (i.e., \code{min(table(tcc$group[, 1])) > 1}) 
                       and \code{"deseq"} when analyzing the count data
                       without replicates 
                       (i.e., \code{min(table(tcc$group[, 1])) == 1})).}
    \item{test.method}{character specifying method for identifying 
                       differentially expressed genes (DEGs) used in STEP2. 
                       Possible values are listed on Details section of 
                       \code{\link{estimateDE}}.
                       The default is \code{"edger"} when
                       analyzing the count data with multiple replicates
                       (i.e., \code{min(table(tcc$group[, 1])) > 1}) 
                       and \code{"deseq"} when analyzing the count data 
                       without replicates 
                       (i.e., \code{min(table(tcc$group[, 1])) == 1))}.}
    \item{iteration}{logical or numeric value specifying the number 
                     of iteration (\eqn{n}) in the proposed normalization 
                     pipeline: the 
                     \eqn{STEP1-(STEP2-STEP3)_{n}}{STEP1-(STEP2-STEP3)n} 
                     pipeline. If \code{FALSE} or \code{0} is specified, 
                     the normalization pipeline is performed only by the 
                     method in STEP1. If \code{TRUE} or \code{1} is
                     specified, the three-step normalization pipeline is 
                     performed. Integers higher than \code{1} indicate the
                     number of iteration in the pipeline.}
    \item{FDR}{numeric value (between 0 and 1) specifying the threshold for
               determining DEGs after STEP2.}
    \item{floorPDEG}{numeric value (between 0 and 1) specifying the 
                     minimum value to be eliminated as potential DEGs 
                     before performing STEP3.}
    \item{floor.value}{data less than \code{floor.value} are placed by 
                       \code{floor.value}. The default is 1.}
    \item{increment}{logical value. if \code{increment = TRUE}, the DEGES
                     pipeline will perform again from the current iterated
                     result.}
    \item{...}{arguments will pass to \code{\link{estimateDE}} for
               identification at STEP2. See \code{\link{estimateDE}} to 
               check the possible arguments.}
}
\details{
The \code{\link{calcNormFactors}} function is the main function
in the TCC package. 
Since this pipeline employs the DEG identification method at STEP2, 
our multi-step strategy can eliminate the biased effect of potential DEGs 
before the second normalization at STEP3. 
To fully utilize the differentially expressed gene elimination 
strategy (DEGES), we strongly recommend not to use 
\code{iteration = 0} or \code{iteration = FALSE}. 
This function internally calls functions implemented in the other R packages
according to the specified parameters. 

\itemize{
    \item \code{norm.method = "tmm"} and \code{test.method = "edger"}\cr
          The \code{\link[edgeR]{calcNormFactors}} function implemented
          in edgeR is used for obtaining the TMM normalization factors
          at both STEP1 and 3.
    \item \code{norm.method = "deseq"}\cr
          The \code{\link[DESeq]{estimateSizeFactors}} function is used
          for obetaining the size factors. The size factors internally 
          convert to normalization factors that are comparable to the 
          TMM normalization factors.
}

The \code{test.method} can be specified \code{"edger"}, \code{"deseq"}, etc.
See DEtails section in \code{\link{estimateDE}}
to check the possible paramaters are used in \code{test.method}.
}
\value{
\link{TCC-class} object containing the normalization factors in the 
\code{norm.factors} field (numerical vector). 
In other words, the normalization factors in the \code{norm.factors} 
field are populated.
The information about normalization, such as execution time, potential
DEGs in STEP2 and DEGES pipeline, are storing in \code{DEGES} field 
(i.e., \code{tcc$DEGES}).
}
\examples{
data(hypoData)
group <- c(1, 1, 1, 2, 2, 2)

# calculating normalization factors using the DEGES/edgeR method 
# (the TMM-edgeR-TMM pipeline)
tcc <- new("TCC", hypoData, group)
tcc <- calcNormFactors(tcc, norm.method = "tmm", test.method = "edger",
                       iteration = 1, FDR = 0.1, floorPDEG = 0.05)
tcc$norm.factors

# calculating normalization factors using the iterative DEGES/edgeR method 
# (iDEGES/edgeR) with n = 3
tcc <- new("TCC", hypoData, group)
tcc <- calcNormFactors(tcc, norm.method = "tmm", test.method = "edger",
                       iteration = 3, FDR = 0.1, floorPDEG = 0.05)
tcc$norm.factors

# calculating normalization factors for simulation data without replicates
tcc <- simulateReadCounts(replicates = c(1, 1))
tcc <- calcNormFactors(tcc, norm.method = "deseq", test.method = "deseq",
                       iteration = 1, FDR = 0.1, floorPDEG = 0.05)
tcc$norm.factors
}
