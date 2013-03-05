\name{calcNormFactors}
\docType{methods}
\alias{calcNormFactors,TCC-method}
\alias{calcNormFactors,DGEList-method}
\alias{calcNormFactors}
\title{Calculate normalization factors}
\description{
This function calculates normalization factors using a specified 
multi-step normalization method from a \link{TCC-class} object.
The procedure can generally be described as the \eqn{STEP1-(STEP2-STEP3)_{n}}{STEP1-(STEP2-STEP3)n} pipeline.
}
\synopsis{
calcNormFactors(tcc, ...)
}
\usage{
calcNormFactors(tcc, norm.method = c("tmm", "deseq"), 
                test.method = c("edger", "deseq", "bayseq"),
                iteration = 1, FDR = 0.1, floorPDEG = 0.05, 
                dispersion =NULL, 
                samplesize = 10000, processors = NULL)
}
\arguments{
	\item{tcc}{\link{TCC-class} object.}
	\item{norm.method}{character specifying normalization method that
	 is used in both the STEP1 and STEP3. Possible values are \code{"tmm"} 
	 for the TMM normalization method implemented in the edgeR package,
	 \code{"edger"} (same as \code{"tmm"}), and \code{"deseq"} for the
	 method implemented in the DESeq package. The default is \code{"tmm"}
	 when analyzing the count data with multiple replicates 
	 (i.e., \code{min(tcc$group) > 1}) and \code{"deseq"} when analyzing
	 the count data with a single measurement (i.e., \code{min(tcc$group) == 1})).}
	\item{test.method}{character specifying method for identifying 
	 differentially expressed genes (DEGs) used in STEP2. Possible
	 values are \code{"edger"}, \code{"deseq"}, and \code{"bayseq"}
	 for the DEG identification methods implemented in the edgeR, DESeq,
	 and baySeq, respectively. The default is \code{"edger"} when
	 analyzing the count data with multiple replicates
	 (i.e., \code{min(tcc$group) > 1}) and \code{"deseq"} when 
	 analyzing the count data with a single measurement (i.e., \code{min(tcc$group) == 1))}.}
	\item{iteration}{logical or numeric value specifying the number of iteration (n) 
	 in the proposed normalization pipeline: 
	 the \eqn{STEP1-(STEP2-STEP3)_{n}}{STEP1-(STEP2-STEP3)n} pipeline.
	 If \code{FALSE} or \code{0} is specified, the normalization pipeline
	 is performed only by the method in STEP1. If \code{TRUE} or \code{1} is
	 specified, the three-step normalization pipeline is performed. Integers
	 higher than \code{1} indicate the number of iteration in the pipeline.}
	\item{FDR}{numeric value (between 0 and 1) specifying the 
	 threshold for determining DEGs after STEP2.}
	\item{floorPDEG}{numeric value (between 0 and 1) specifying the minimum value 
	 to be eliminated as potential DEGs before performing STEP3.}
	\item{dispersion}{either a numeric vector of dispersions indicating that
	dispersions. See the \code{\link[edgeR]{exactTest}} function for details.}
	\item{samplesize}{numeric value specifying the sample size for estimating
	 the prior parameters if \code{test.method = "bayseq"}. See 
	 the \code{\link[baySeq]{getPriors.NB}} function for details.}
	\item{processors}{numeric value or 'snow' object for using multi 
	 processors if \code{test.method = "bayseq"}. See the 
	\code{\link[baySeq]{getPriors.NB}} function for details.}
}
\details{
The \code{\link{calcNormFactors}} function is the main function in the TCC package. 
Since this pipeline employs the DEG identification method at STEP2, 
our multi-step strategy can eliminate the biased effect of potential DEGs 
before the second normalization at STEP3. 
To fully utilize the differentially expressed gene elimination strategy (DEGES), 
we strongly recommend not to use \code{iteration = 0} or \code{iteration = FALSE}. 
This function internally calls functions implemented in the edgeR, 
DESeq, and baySeq packages according to the specified parameters. 

If the \code{norm.method = "tmm"} is specified, 
the \code{calcNormFactors} function implemented in edgeR 
is used for obtaining the TMM normalization factors at both STEP1 and 3. 
In case of \code{norm.method = "deseq"}, 
the \code{\link[DESeq]{estimateSizeFactors}} function in DESeq is used. 
Note that the original \code{\link[DESeq]{estimateSizeFactors}} function 
returns the size factors (not normalization factors). 
Our \code{\link{calcNormFactors}} function internally converts the size factors 
into normalization factors that are comparable to the TMM normalization factors.

If the \code{test.method = "edger"} is specified, 
a series of functions for differential expression 
analysis (\code{\link[edgeR]{estimateCommonDisp}}, 
\code{\link[edgeR]{estimateTagwiseDisp}}, 
and \code{\link[edgeR]{exactTest}}) in edgeR are internally used. 
Similarly, the \code{test.method = "deseq"} internally use two functions 
(\code{\link[DESeq]{estimateDispersions}} and \code{\link[DESeq]{nbinomTest}}) 
in DESeq and the \code{test.method = "bayseq"}
 internally use two functions (\code{\link[baySeq]{getPriors.NB}} 
and \code{\link[baySeq]{getLikelihoods.NB}}) in baySeq.
}
\value{
\link{TCC-class} object containing the normalization factors in the 
\code{norm.factors} field (numerical vector). 
In other words, the normalization factors in the \code{norm.factors} field are populated.
}
\examples{
# calculating normalization factors for a hypothetical count data 
# using the TbT method (the TMM-baySeq-TMM pipeline; Kadota et al., 2012)
\dontrun{
data(hypoData)
group <- c(3, 3)
tcc <- new("TCC", hypoData, group)
tcc <- calcNormFactors(tcc, norm.method = "tmm", test.method = "bayseq")
tcc$norm.factors
}

# calculating normalization factors using the DEGES/edgeR method 
# (the TMM-edgeR-TMM pipeline)
\dontrun{
tcc <- new("TCC", hypoData, group)
tcc <- calcNormFactors(tcc)
tcc$norm.factors
}

# calculating normalization factors using the iterative DEGES/edgeR method 
# (iDEGES/edgeR) with n = 3
\dontrun{
tcc <- new("TCC", hypoData, group)
tcc <- calcNormFactors(tcc, iteration = 3)
tcc$norm.factors
}

# calculating normalization factors for count data without replicates
\dontrun{
group <- c(1, 1)
tcc <- new("TCC", hypoData[, c(1, 4)], group)
tcc <- calcNormFactors(tcc)
tcc$norm.factors
}

# calculating normalization factors for another simulation data without 
# replicates
\dontrun{
tcc <- generateSimulationData(group = c(1, 1))
tcc <- calcNormFactors(tcc)
tcc$norm.factors
}

# calculating normalization factors for a simulation data using the 
# TbT method with four processors
\dontrun{
tcc <- generateSimulationData()
tcc <- calcNormFactors(tcc, norm.method = "tmm", 
                       test.method = "bayseq", processors = 4)
tcc$norm.factors
}

# calculating normalization factors for a simulation data using 
# the TbT method with four processors (another way)
\dontrun{
tcc <- generateSimulationData()
cl <- makeCluster(2, "SOCK")
tcc <- calcNormFactors(tcc, norm.method = "tmm", 
                       test.method = "bayseq", processors = cl)
tcc$norm.factors
}
}
