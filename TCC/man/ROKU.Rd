\name{ROKU}
\alias{ROKU}
\title{detect tissue-specific genes in microarray data}
\description{
This function detects tissue-specific genes from
a matrix or data frame object.
}
\usage{
ROKU(data, upper.limit = 0.25, sort = FALSE)
}
\arguments{
    \item{data}{numeric matrix containing microarray data, where row
                indicates the gene or probeset ID, each column indicates
                the tissue, and each cell indicates the data of the gene
                in the tissue.}
    \item{upper.limit}{the percentage of all tissues specifying
                       the upper limit for detecting outliers.}
    \item{sort}{logical. If \code{TRUE}, the retrieved results are sorted
                in order of the rank of entropy score. If \code{FALSE}, the
                results are retrieved by the original order.}
}
\details{
ROKU is a method for detecting tissue-specific genes from microarry data.
In detailed, ROKU analyzes microarray data in two steps.
The first step is to compute the Shannon entropy for all genes,
and ranks genes according to overall tissue-specificity with it.
The second step is to identify the specific tissues whose observations
are regarded as outliers using Akaike's information criterion procedure.

The \code{ROKU} function retrieves a data frame containing
outlier matrix, the entropy score, and the rank of entropy score.
}
\value{
A data frame consisting an outiler matrix,
entropy scores, and the rank of entropy scores.

The last two columns of retrived results indicate 
the entropy scores and it's ranks.
The columns from the first to the third-to-last indicate
the outlier matrix consisting of 1, -1, and 0.
1 for over-expressed outliers,
-1 for under-expressed outliers,
and 0 for non-expressed outliers.
}
\references{
Kadota K, Ye J, Nakai Y, Terada T, Shimizu K: ROKU: a novel method for identification of tissue-specific genes. BMC Bioinformatics 2006, 7: 294.

Kadota K, Nishimura SI, Bono H, Nakamura S, Hayashizaki Y, Okazaki Y, Takahashi K: Detection of genes with tissue-specific expression patterns using Akaike's Information Criterion (AIC) procedure. Physiol Genomics 2003, 12:251-259.

Ueda T. Simple method for the detection of outliers. Japanese J Appl Stat 1996, 25: 17-26.
}
\examples{
# create sample data.
x <- matrix(rnorm(100), ncol = 10)

roku <- ROKU(x)
}
