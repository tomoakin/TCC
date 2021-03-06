\name{hypoData_mg}
\docType{data}
\alias{hypoData_mg}
\title{A simulation dataset for comparing two-group tag count
 data, focusing on RNA-seq}
\description{
A simulation dataset, consisting of 1,000 rows, corresponding to genes, and
 9 columns, representing 3 independent biological samples from 3 conditions.
}
\usage{data(hypoData_mg)}
\format{hypoData is a matrix with 9 columns and 1,000 rows.}
\details{
This package typically start the differential expression analysis with
 a count table matrix such as hypoData_mg where each row indicates the
 gene (or transcript), each column indicates the sample (or library),
 and each cell indicates the number of counts to the gene in the sample.
 The first three columns are produced from biological replicates of,
 for example, Group 1, the next three columns are from Group2 and the
 remaining columns are from Group 3;
 i.e., {G1_rep1, G1_rep2, G1_rep3} vs. {G2_rep1, G2_rep2, G2_rep3} vs.
 {G3_rep1, G3_rep2, G3_rep3}.
 This data is generated by the \code{\link{simulateReadCounts}}
 function with default parameter settings. The first 200 genes are
 differentially expressed among the three groups. Of these, the first 140
 genes are expressed at a higher level only in Group 1 (G1), the next
 40 genes are expressed at a higher level only in G2 and the remaining
 20 genes are expressed at a higher level only in G3. Accordingly, the
 201-1000th genes are not differentially expressed (non-DEGs).
 The levels of differential expression (DE) are four-fold in both groups. 

# The 'hypoData_mg' is generated by following commands.
tcc <- simulateReadCounts(Ngene = 1000, PDEG = 0.2,
                              DEG.assign = c(0.7, 0.2, 0.1),
                              DEG.model = "uniform",
                              DEG.foldchange = c(4, 4),
                              replicates = c(3, 3, 3))
hypoData_mg <- tcc$count
}
\keyword{datasets}
