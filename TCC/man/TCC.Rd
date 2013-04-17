\name{TCC}
\docType{package}
\alias{TCC}
\alias{TCC-package}
\title{A package for differential expression analysis from two-group tag
 count data with robust normalization strategies}
\description{
This package performs differential expression analysis from tag count 
data that are produced from high-throughput sequencing (HTS) or next 
generation sequencing (NGS) technology. A notable feature of this 
package is to provide robust normalization methods whose strategy 
is to remove data assigned as potential differentially expressed 
genes (DEGs) before performing data normalization (Kadota et al., 2012).
}
\details{
TCC is a package for differential expression analysis from tag count data,
focusing of RNA-seq. This package implements some functions 
for calculating normalization factors,
identifying DEGs,
depicting so-called M-A plot,
and generating simulation data.

To utilize this package, 
the count matrix coupled with label information should be stored to a
 \link{TCC-class} object using the \code{new} method.
All functions (except for some legacy functions) used in this package
require this \link{TCC-class} object.
Using this object, the \code{\link{calcNormFactors}} function calculates
normalization factors and the \code{\link{estimateDE}} function estimates
the degree of differential expression (DE) for individual genes. 
The estimated normalization factors obtained by using the
\code{\link{calcNormFactors}} function are used within the statistical
model for differential analysis in the \code{\link{estimateDE}} function.
Both two functions internally call functions from other packages
(edgeR, DESeq, and baySeq) when needed. TCC also provides a useful function
\code{\link{simulateReadCounts}} for generating simulation data 
with various conditions.
Convenient plotting functions are also included.
}
\seealso{
\link{TCC-class}
}
\references{
Robinson MD, McCarthy DJ, Smyth GK: edgeR: a Bioconductor package for
 differential expression analysis of digital gene expression data.
 Bioinformatics 2010, 26(1): 139-140

Hardcastle TJ, Kelly KA: baySeq: empirical Bayesian methods for
 identifying differential expression in sequence count data.
 BMC Bioinformatics 2010, 11: 422

Anders S, Huber W: Differential expression analysis for sequence count data.
 Genome Biol. 2010, 11: R106

Kadota K, Nishiyama T, Shimizu K: A normalization strategy for comparing
 tag count data. Algorithms Mol Biol. 2012, 7:5
}
\keyword{packages}
