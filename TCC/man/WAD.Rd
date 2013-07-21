\name{WAD}
\alias{WAD}
\title{Calculate WAD statistic for individual genes}
\description{
This function caclulates WAD statistic from a matrix or
data.frame object.
}
\usage{
WAD(data, group, log.scale = FALSE, floor.value = 1, sort = FALSE)
}
\arguments{
    \item{data}{numeric matrix containing count data or 
                microarray data. Each row indicates the gene or probeset ID,
                each column indicates the sample or group, and each cell
                indicates the data of the gene in the sample.}
    \item{group}{numeric vector indicating the experimental group for each
                 sample (or group).}
    \item{log.scale}{data using in \code{WAD} function should be log scaled.
                     If \code{TRUE}, data convert to log scale with base of 2.
                     The default is \code{FALSE}.}
    \item{floor.value}{if log scaled data less than \code{floor.value} are
                       placed by \code{floor.value}. The default is 1. Ignored
                       if \code{log.scale} is \code{FALSE}.}
    \item{sort}{logical. If \code{TRUE}, the retrieved results are sorted
                in order of the rank of absolute WAD statistic.
                If \code{FALSE}, the results are retrieved by the original
                order.}
}
\details{
WAD is a method for ranking differentially expressed genes (DEGs)
for two-group data from a \code{matrix} or \code{data.frame} object.
The method reports WAD statistic for individual genes.
The large absolute value of WAD statsitic is a evidence
of a differential expression.

\code{WAD} function calculates the WAD statistic from a \code{matrix} 
or \code{data.frame} obeject. The data used in \code{WAD} function
need convert to log scale with base of 2.
The \code{WAD} function provides \code{log.scale} argument for 
automatically log scaling.
}
\value{
A numeric vector of WAD statistic for individual genes
}
\references{
Kadota K, Nakai Y, Shimizu K: A weighted average difference method for
detecting differentially expressed genes from microarray data.
Algorithms Mol Biol. 2008, 3:8.
}
\examples{
data(hypoData)
data <- log2(hypoData)
group <- c(1, 1, 1, 2, 2, 2)

wad <- WAD(data, group, sort = TRUE)
}
