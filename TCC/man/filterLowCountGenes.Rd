\name{filterLowCountGenes}
\alias{filterLowCountGenes}
\title{Filter genes from a TCC-class object}
\description{
This function takes a TCC object and returns a new TCC object
without genes having low count tags across samples. 
The threshold is configurable with \code{low.count} parameter.
}
\usage{filterLowCountGenes(tcc, low.count = 0)}
\arguments{
	\item{tcc}{\link{TCC-class} object.}
	\item{low.count}{numeric value (>= 0) specifying the threshold
	 for filtering genes. The higher value indicates the more
	 numbers of genes to be filtered out.}
}
\value{
\link{TCC-class} object consisting of genes whose sum of the counts
	 across samples is equal or higher than the \code{low.count} value.
}
\examples{
# Filtering genes with zero counts across samples (default) from
# a hypothetical count dataset that originally has 10,000 genes
data(hypoData)
group <- c(1, 1, 1, 2, 2, 2)
tcc <- new("TCC", hypoData, group)
dim(tcc$count)
tcc0 <- filterLowCountGenes(tcc)
dim(tcc0$count)

# Filtering genes with 10 counts across samples from hypoData
dim(tcc$count)
tcc1 <- filterLowCountGenes(tcc, low.count = 10)
dim(tcc1$count)
}
\keyword{methods}
