\name{plotFCPseudocolor}
\alias{plotFCPseudocolor}
\title{Create a pseudo-color image of simulation data}
\description{
This function creates a pseudo-color image of simulation data 
regarding the number of differentially expressed genes (DEGs) 
and the breakdowns for individual groups from a \link{TCC-class} object.
}
\usage{
plotFCPseudocolor(tcc, main, xlab, ylab)
}
\arguments{
	\item{tcc}{\link{TCC-class} object.}
	\item{main}{character string indicating the plotting title.}
	\item{xlab}{character string indicating the \emph{x}-label title.}
	\item{ylab}{character string indicating the \emph{y}-label title.}
}
\details{
This function should be used after the 
\code{\link{generateSimulationData}} function that generates 
simulation data with arbitrary defined conditions. 
The largest log fold-change (FC) values are 
in magenta and no-changes are in white. 
}
\examples{
# Generating a simulation data for comparing two groups
# (G1 vs. G2) with biological replicates.
# the first 20 genes are DEGs, where 18 are up in G1.
tcc <- generateSimulationData(Ngene = 100, PDEG = 0.2,
                         DEG.assign = c(0.9, 0.1),
                         DEG.model = "uniform",
                         DEG.foldchange = list(4, 4), replicates = c(3, 3))
plotFCPseudocolor(tcc)

# Generating a simulation data for comparing three groups
# (G1 vs. G2 vs. G3) with biological replicates.
# the first 3000 genes are DEGs, where the 70%, 20%, and 10% are
# up-regulated in G1, G2, G3, respectively. The levels of DE are
# 3-, 10, and 6-fold in individual groups.
tcc <- generateSimulationData(Ngene = 10000, PDEG = 0.3,
                         DEG.assign = c(0.7, 0.2, 0.1),
                         DEG.model = "uniform",
                         DEG.foldchange = list(3, 10, 6),
                         replicates = c(2, 3, 3))
plotFCPseudocolor(tcc)
}
\keyword{methods}
