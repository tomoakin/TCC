\name{nakai}
\docType{data}
\alias{nakai}
\title{DNA microarray data set}
\description{
The DNA microarray data 
}
\usage{data(nakai)}
\format{\code{nakai} is a matrix of dimension 32,099 times 8.}
\details{
This is a DNA microarray data set which consists of 32,099 probes
of eight samples.
The eight samples are taken from the liver which were obtained from 
rats that had been fed and fasted for 24 hours.

The cDNA was synthesized from purified total RNA,
and added to an Affymetrix Rat Genome 230 2.0 Array.
The fluorescence signals were scanned by using
Affymetrix GeneChip System.

The signals saved in \code{nakai} data set are rounded
of to three decimal places. 
The original data set are available on NCBI Gene Expression Omnibus
with GSE7623.
} 
\references{
Nakai Y, Hashida H., Kadota K, Minami M, Shimizu K, Matsumoto I, Kato H, Abe K., Up-regulation of genes related to the ubiquitin-proteasome system in the brown adipose tissue of 24-h-fasted rats. Bioscience Biotechnology and Biochemistry 2008, 72(1):139-148.
}
\examples{
data(nakai)
}
\keyword{datasets}
