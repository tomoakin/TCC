\name{estimateDE}
\alias{estimateDE}
\title{Estimate degrees of differential expression (DE) for individual genes}
\usage{
estimateDE(tcc, test.method = NULL, FDR = NULL, 
           dispersion = NULL, fit0 = NULL, fit1 = NULL, design = NULL,
           contrast = NULL, coef = NULL, comparison = NULL,
           samplesize = NULL, floor.value = 1, cl = NULL)
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
                       DEGs. Possible values are listed on Details section.
                       The default is \code{"edger"} when analyzing the 
                       count data with replicates (i.e.,
                       \code{min(table(tcc$group[, 1])) > 1}) and 
                       \code{"deseq"} when analyzing the 
                       count data without replicates 
                       (i.e., \code{min(table(tcc$group[, 1])) == 1)}).}
    \item{FDR}{numeric value (between 0 and 1) specifying the threshold
               for determining DEGs.}
    \item{dispersion}{numeric vector giving the dispersion of all genes 
                      for analyzing the count data without replicates when 
                      the \code{test.method = "edger"} is specified.
                      See the \code{\link[edgeR]{exactTest}} function in 
                      edgeR for details.}
    \item{design}{numeric matrix giving the design matrix for the linear 
                  model. See the \code{\link[edgeR]{glmFit}} function in
                  edgeR for details.}
    \item{contrast}{numeric vector specifying a contrast of the linear
                    model coefficients to be tested equal to zero. See the 
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
    \item{samplesize}{intenger specifying the sample size for
                       estimating the prior parameters if 
                      \code{test.method = "bayseq"} (default is 10000).
                      If \code{test.method = "ebseq"}, the value is used for
                      specifying the ietration of EM algorithm in EBSeq 
                      (default is 5).
                      If \code{test.method = "samseq"}, the value is used for
                      specifying the number of permulations in samr 
                      (default is 10).}
    \item{floor.value}{data less than \code{floor.value} are placed by 
                       \code{floor.value}. The defaul is 1.}
    \item{cl}{'snow' object for using multi processors if 
              \code{test.method = "bayseq"}. See the 
              \code{\link[baySeq]{getPriors.NB}} function for details.}
}
\details{
\code{estimateDE} function is generally used after the 
\code{\link{calcNormFactors}} function calculated normalization factors. 
\code{estimateDE} constructs a statistical model for differential 
expression (DE) analysis with the calculated normalization factors. 
This function internally calls individual functions implemented in the 
other R packages according to the specified parameters.

\itemize{
    \item \code{test.method = "edger"}\cr
          \code{\link[edgeR]{estimateCommonDsip}},
          \code{\link[edgeR]{estimateTagwiseDisp}},
          \code{\link[edgeR]{exactTest}} in edgeR are 
          internally called for two-group comparison. \cr
          \code{\link[edgeR]{estimateGLMCommonDisp}},
          \code{\link[edgeR]{estimateGLMTrendedDisp}},
          \code{\link[edgeR]{estimateGLMTagwiseDisp}},
          \code{\link[edgeR]{glmFit}},
          \code{\link[edgeR]{glmLRT}}
          are for multi-group and multi-factor comparison.
    \item \code{test.method = "deseq"}\cr
          \code{\link[DESeq]{estimateDispersions}},
          \code{\link[DESeq]{nbinomTest}} in DESeq are
          internally called for two-group comparison. \cr
          \code{\link[DESeq]{estimateDispersions}},
          \code{\link[DESeq]{fitNbinomGLMs}},
          \code{\link[DESeq]{nbinomGLMTest}}
          are for multi-group and multi-factor comparison.
    \item \code{test.method = "bayseq"}\cr
          \code{\link[baySeq]{getPriors.NB}},
          \code{\link[baySeq]{getLikelihoods.NB}}) in baySeq
          are called.
          Since baySeq outputs posterior likelihoods instead of
          \eqn{p}-values, the \eqn{(1 - likelihood)} values are
          stored in the corersponding field
          (i.e., \code{tcc$stat$p.value}).
    \item \code{test.method = "noiseq"}\cr
          \code{\link[NOISeq]{noiseq}} in NOISeq are called with
          the paramaters of \code{k = 0.5} and
          \code{repliactes = "biological"}.
          Since NOISeq outputs probability of DEGs, the probability are
          stored in \code{tcc$stat$prob} field, 
          and \code{tcc$stat$p.value} be \code{NULL}.
    \item \code{test.method = "ebseq"}\cr
          \code{\link[EBSeq]{EBTest}} in EBSeq is called 
          for two-group comparison. \cr
          \code{\link[EBSeq]{EBMultiTest}} is called for 
          multi-group comparison.
          Since EBSeq only outputs \eqn{q}-value,
          the all values in \code{tcc$stat$p.value} are \code{NA}.
    \item \code{test.method = "samseq"}\cr
          \code{\link[samr]{SAMseq}} in samr is called.
          Since \code{\link[samr]{SAMseq}} function outputs test statistics
          instead of \eqn{p}-values, 
          the \code{tcc$stat$p.value} be \code{NULL}.
          The test statistics are stored in \code{tcc$stat$testStat} filed.
    \item \code{test.method = "wad"}
          The \code{\link{WAD}} function is used.
          Since \code{\link[samr]{SAMseq}} function outputs test statistics
          instead of \eqn{p}-values, 
          the \code{tcc$stat$p.value} be \code{NULL}.
          The test statistics are stored in \code{tcc$stat$testStat} filed.
}
}
\value{
A \code{\link{TCC-class}} object containing following fields:
    \item{stat$p.value}{numeric vector of \eqn{p}-values.}
    \item{stat$q.value}{numeric vector of \eqn{q}-values calculated
                        based on the \eqn{p}-values using the 
                        \code{p.adjust} function
                        with default parameter settings.}
     \item{stat$rank}{gene rank in order of the \eqn{p}-values.}
     \item{estimatedDEG}{numeric vector consisting of 0 or 1
                         depending on whether each gene is classified 
                         as non-DEG or DEG. The threshold for classifying
                         DEGs or non-DEGs is preliminarily given as the 
                         \code{FDR} argument.}
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
tcc <- calcNormFactors(tcc, norm.method = "tmm", test.method = "edger",
                       iteration = 1, FDR = 0.1, floorPDEG = 0.05)
tcc <- estimateDE(tcc, test.method = "edger", FDR = 0.1)
head(tcc$stat$p.value)
head(tcc$stat$q.value)
head(tcc$estimatedDEG)
result <- getResult(tcc)


# Analyzing a simulation data for comparing two groups
# (G1 vs. G2) without replicates
# The DE analysis is performed by an negative binomial (NB) test
# in DESeq coupled with the DEGES/DESeq normalization factors.
data(hypoData)
group <- c(1, 2)
tcc <- new("TCC", hypoData[, c(1, 4)], group)
tcc <- calcNormFactors(tcc, norm.method = "deseq", test.method = "deseq",
                       iteration = 1, FDR = 0.1, floorPDEG = 0.05)
tcc <- estimateDE(tcc, test.method = "deseq", FDR = 0.1)
}
\keyword{methods}
