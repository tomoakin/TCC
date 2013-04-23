\name{estimateDE}
\alias{estimateDE}
\title{Estimate degrees of differential expression (DE) for individual genes}
\usage{
estimateDE(tcc, test.method = NULL, FDR = NULL,  
           samplesize = 10000, cl = NULL, \dots)
}
\description{
This function calculates \eqn{p}-values (or the related statistics) for
identifying differentially expressed genes (DEGs) from a
\link{TCC-class} object. 
\code{estimateDE} internally calls a specified method 
implemented in other R packages.
}
\arguments{
	\item{tcc}{\link{TCC-class} object.}
	\item{test.method}{character string specifying method for identifying
	 DEGs. Possible values are \code{"edger"}, \code{"deseq"}, and
	 \code{"bayseq"} for the DEG identification methods implemented in
	 the edgeR, DESeq, and baySeq, respectively. The default is
	 \code{"edger"} when analyzing the count data with replicates (i.e.,
	 \code{min(tcc$group) > 1}) and \code{"deseq"} when analyzing the 
	 count data without any replicates (i.e., \code{min(tcc$group) == 1)}).}
	\item{FDR}{numeric value (between 0 and 1) specifying the threshold
	 for determining DEGs.}
	\item{samplesize}{numeric value specifying the sample size for
	 estimating the prior parameters if \code{test.method = "bayseq"}.
	 See the \code{\link[baySeq]{getPriors.NB}} function for details.}
	\item{cl}{cluster object of 'snow' package for using multiple
	 processors. Effective when \code{test.method = "bayseq"}. See the
	 \code{\link[baySeq]{getPriors.NB}} function for details.}
        \item{\dots}{parameters to pass to edgeR, DESeq, or baySeq.}
}
\details{
\code{estimateDE} function is generally used after the 
\code{\link{calcNormFactors}} function calculated normalization factors. 
\code{estimateDE} constructs a statistical model for differential 
expression (DE) analysis with the calculated normalization factors. 
This function internally calls individual functions implemented in the edgeR, 
DESeq, and baySeq packages according to the specified parameters.

If the \code{test.method = "edger"} is specified, 
a series of functions for differential expression analysis 
(\code{\link[edgeR]{estimateCommonDisp}},
 \code{\link[edgeR]{estimateTagwiseDisp}}, 
and \code{\link[edgeR]{exactTest}}) 
in edgeR are internally called and \eqn{p}-values (and the derivative such
 as the \eqn{q}-values and the ranks) are calculated. 

The \code{test.method = "deseq"} internally use two functions 
(\code{\link[DESeq]{estimateDispersions}} and \code{\link[DESeq]{nbinomTest}})
 in DESeq. 

The \code{test.method = "bayseq"} internally use two 
functions (\code{\link[baySeq]{getPriors.NB}} and
 \code{\link[baySeq]{getLikelihoods.NB}}) in baySeq.

Different from the edgeR and DESeq in which the calculated \eqn{p}-values are
stored in the \code{stat$p.value} field of the \link{TCC-class} object
\code{tcc}, baySeq outputs posterior likelihoods instead of \eqn{p}-values. 
Therefore, the \eqn{(1 - likelihood)} values are stored in the corresponding
 field in case of \code{test.method = "bayseq"}.
}
\value{
A \code{\link{TCC-class}} object containing following fields:
	\item{stat$p.value}{numeric vector of \eqn{p}-values.}
	\item{stat$q.value}{numeric vector of \eqn{q}-values calculated
	 based on the \eqn{p}-values using the \code{p.adjust} function
	 with default parameter settings.}
	\item{stat$rank}{gene rank in order of the \eqn{p}-values.}
	\item{estimatedDEG}{numeric vector consisting of 0, 1, or 2
	 depending on whether each gene is classified as non-DEG, DEG
	 expressed at a higher level in Group 1 or Group 2, respectively.
	 The threshold for classifying DEGs or non-DEGs is preliminarily
	 given as the \code{FDR} argument.}
}
\examples{
# Analyzing a simulation data for comparing two groups
# (G1 vs. G2) with biological replicates
# The DE analysis is performed by an exact test in edgeR coupled
# with the DEGES/edgeR normalization factors.
# For retrieving the summaries of DE results, we recommend to use
# the getResult function.
data(hypoData)
hypoData <- hypoData[1:100*10,]
group <- c(1, 1, 1, 3, 3, 3)
tcc <- new("TCC", hypoData, group)
tcc <- calcNormFactors(tcc)
tcc <- estimateDE(tcc, test.method = "edger", FDR = 0.1)
head(tcc$stat$p.value)
head(tcc$stat$q.value)
head(tcc$estimatedDEG)


# Analyzing a simulation data for comparing two groups
# (G1 vs. G2) without any replicates
# The DE analysis is performed by an negative binomial (NB) test
# in DESeq coupled with the DEGES/DESeq normalization factors.
group <- c(1, 2)
tcc <- new("TCC", hypoData[, c(1, 4)], group)
tcc <- calcNormFactors(tcc)
tcc <- estimateDE(tcc, test.method = "deseq", FDR = 0.1)
}
\keyword{methods}
