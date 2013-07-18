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
\S4method{calcNormFactors}{TCC}(tcc, norm.method = NULL, 
                test.method = NULL, iteration = TRUE, 
                FDR = NULL, floorPDEG = 0.05, 
                dispersion = NULL, 
                design = NULL, contrast = NULL, coef = NULL,
                fit0 = NULL, fit1 = NULL, comparison = NULL,
                samplesize = 10000, cl = NULL, 
                trimWAD = 0.50, q = 0.25, AD = FALSE, k = 1,
                increment = FALSE)
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
                       Possible values are \code{"edger"}, \code{"deseq"}, 
                       and \code{"bayseq"} for the DEG identification methods
                       implemented in the edgeR, DESeq, and baySeq, 
                       respectively. The default is \code{"edger"} when
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
    \item{dispersion}{numeric vector giving the dispersion of all genes for
                      analyzing the count data without replicates when the 
                      \code{test.method = "edger"} is specified.
                      See the \code{\link[edgeR]{exactTest}} function in 
                      edgeR for details.}
    \item{design}{numeric matrix giving the design matrix for the linear 
                  model. See the \code{\link[edgeR]{glmFit}} function in 
                  edgeR for details.}
    \item{contrast}{numeric vector specifying a contrast of the linear model
                    coefficients to be tested equal to zero. See the 
                    \code{\link[edgeR]{glmLRT}} function in edgeR for 
                    details.}
    \item{coef}{integer or character vector indicating which coefficients
                of the linear model are to be tested equal to zero. See the
                \code{\link[edgeR]{glmLRT}} function in edgeR for details.}
    \item{fit0}{a formula for creating reduced model described in DESeq. 
                The left hand side must be 'count', the right hand side 
                can involve any column of \code{tcc$group} is used as the 
                model frame. See the \code{\link[DESeq]{fitNbinomGLMs}} 
                function for details.}
    \item{fit1}{a formula for creating full model described in DESeq.
                The left hand side must be 'count', the right hand side 
                can involve any column of \code{tcc$group} is used as the 
                model frame. See the \code{\link[DESeq]{fitNbinomGLMs}} 
                function for details.}
    \item{comparison}{numeric or character string identifying the columns 
                      in the \code{tcc$group[, 1]} for analysis. 
                      See the \code{group} argument in 
                      \code{\link[baySeq]{topCounts}} for details.}
    \item{samplesize}{numeric value specifying the sample size for estimating
                      the prior parameters if \code{test.method = "bayseq"}. 
                      See the \code{\link[baySeq]{getPriors.NB}} function 
                      for details.}
    \item{k}{data equal to zero are placed by k.}
    \item{cl}{'snow' object for using multi processors if 
              \code{test.method = "bayseq"}. See the 
              \code{\link[baySeq]{getPriors.NB}} function for details.}
    \item{q}{numeric value (between 0 and 1) specifying the 
             \eqn{q \times 100}{q X 100}th percentile of the count data
             as a threshold. The counts less than the threshold are trimmed
             before perform TWAD normalization.}
    \item{trimWAD}{numeric value (between 0 and 1) specifying the rate of
                   genes to be trimmed by WAD statistic before perform 
                   TWAD normalization.}
    \item{AD}{logical value specifying the algorithm for calculating WAD
              statistics.}
    \item{increment}{logical value. if \code{increment = TRUE}, the DEGES
                     pipeline will perform again from the current iterated
                     result. For example, assume that the iDEGES pipeline
                     (\code{iteration = n}) was performed, the pipeline
                     performs again if \code{increment = TRUE}. The result
                     is same to (\code{iteration = n + 1}).}
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
This function internally calls functions implemented in the edgeR, 
DESeq, and baySeq packages according to the specified parameters. 

If the \code{norm.method = "tmm"} is specified, 
the \code{calcNormFactors} function implemented in edgeR 
is used for obtaining the TMM normalization factors at both STEP1 and 3. 
In case of \code{norm.method = "deseq"}, 
the \code{\link[DESeq]{estimateSizeFactors}} function in DESeq is used. 
Note that the original \code{\link[DESeq]{estimateSizeFactors}} function 
returns the size factors (not normalization factors). 
Our \code{\link{calcNormFactors}} function internally converts 
the size factors into normalization factors that are comparable 
to the TMM normalization factors.

If the \code{test.method = "edger"} is specified, 
a series of functions for differential expression 
analysis (\code{\link[edgeR]{estimateCommonDisp}}, 
\code{\link[edgeR]{estimateTagwiseDisp}}, 
\code{\link[edgeR]{exactTest}}, etc.) in edgeR are internally used. 
Similarly, the \code{test.method = "deseq"} internally use several functions 
(\code{\link[DESeq]{estimateDispersions}}, 
\code{\link[DESeq]{nbinomTest}}, etc.)
in DESeq and the \code{test.method = "bayseq"}
 internally use two functions (\code{\link[baySeq]{getPriors.NB}} 
and \code{\link[baySeq]{getLikelihoods.NB}}) in baySeq.
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
