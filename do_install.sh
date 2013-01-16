#!/bin/sh
R --slave <<EOI
source("http://bioconductor.org/biocLite.R")
biocLite("BiocUpgrade")
biocLite(c("edgeR", "baySeq", "DESeq", "ROC"))
EOI
time R CMD INSTALL TCC
