\name{estimateDE}
\alias{estimateDE}
\title{Estimate degrees of differential expression (DE) for individual genes}
\usage{
estimateDE(tcc, test.method = NULL, FDR = NULL, 
           dispersion = NULL, fit0 = NULL, fit1 = NULL, design = NULL,
           contrast = NULL, coef = NULL, comparison = NULL,
           samplesize = 10000, cl = NULL)
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
	 \code{min(table(tcc$group[, 1])) > 1}) and \code{"deseq"} when analyzing the 
	 count data without replicates (i.e., \code{min(table(tcc$group[, 1])) == 1)}).}
	\item{FDR}{numeric value (between 0 and 1) specifying the threshold
	 for determining DEGs.}
	\item{dispersion}{numeric vector giving the dispersion of all genes for analyzing
	 the count data without replicates when the \code{test.method = "edger"} is specified.
	 See the \code{\link[edgeR]{exactTest}} function in edgeR for details.}
	\item{design}{numeric matrix giving the design matrix for the linear model.
	 See the \code{\link[edgeR]{glmFit}} function in edgeR for details.}
	\item{contrast}{numeric vector specifying a contrast of the linear
	model coefficients to be tested equal to zero. See the 
	\code{\link[edgeR]{glmLRT}} function in edgeR for details.}
	\item{coef}{integer or character vector indicating which coefficients
	of the linear model are to be tested equal to zero. See the
	\code{\link[edgeR]{glmLRT}} function in edgeR for details.}
	\item{fit0}{a formula for creating reduced model described in DESeq. 
	The left hand side must be 'count', the right hand side can involve 
	any column of \code{tcc$group} is used as the model frame. 
	See the \code{\link[DESeq]{fitNbinomGLMs}} function for details.}
	\item{fit1}{a formula for creating full model described in DESeq.
	The left hand side must be 'count', the right hand side can involve 
	any column of \code{tcc$group} is used as the model frame. 
	See the \code{\link[DESeq]{fitNbinomGLMs}} function for details.}
	\item{comparison}{numeric or character string identifying the columns 
	in the \code{tcc$group[, 1]} for analysis. See the \code{group} argument in 
	\code{\link[baySeq]{topCounts}} for details.}
	\item{samplesize}{numeric value specifying the sample size for
	 estimating the prior parameters if \code{test.method = "bayseq"}.
	 See the \code{\link[baySeq]{getPriors.NB}} function for details.}
	\item{cl}{'snow' object for using multi
	 processors if \code{test.method = "bayseq"}. See the
	 \code{\link[baySeq]{getPriors.NB}} function for details.}
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
 \code{\link[edgeR]{exactTest}}, etc.) 
in edgeR are internally called and \eqn{p}-values (and the derivative such
 as the \eqn{q}-values and the ranks) are calculated. 

The \code{test.method = "deseq"} internally use several functions 
(\code{\link[DESeq]{estimateDispersions}}, \code{\link[DESeq]{nbinomTest}}, etc.)
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
	\item{estimatedDEG}{numeric vector consisting of 0 or 1
	 depending on whether each gene is classified as non-DEG or DEG.
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
group <- c(1, 1, 1, 2, 2, 2)
tcc <- new("TCC", hypoData, group)
tcc <- calcNormFactors(tcc)
tcc <- estimateDE(tcc, test.method = "edger", FDR = 0.1)
head(tcc$stat$p.value)
head(tcc$stat$q.value)
head(tcc$estimatedDEG)


# Analyzing a simulation data for comparing two groups
# (G1 vs. G2) without replicates
# The DE analysis is performed by an negative binomial (NB) test
# in DESeq coupled with the DEGES/DESeq normalization factors.
data(hypoData)
group <- c(1, 2)
tcc <- new("TCC", hypoData[, c(1, 4)], group)
tcc <- calcNormFactors(tcc)
tcc <- estimateDE(tcc, test.method = "deseq", FDR = 0.1)
}
\keyword{methods}
