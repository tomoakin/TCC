set.seed(1000)
library(TCC)
tcc <- simulateReadCounts(Ngene = 20000, PDEG = 0.30,
                              DEG.assign = c(0.85, 0.15),
                              DEG.foldchange = list(min = c(1.2, 1.2),
                                                    shape = c(2.0, 2.0),
                                                    scale = c(0.5, 0.5)),
                              replicates = c(2, 2))
tcc <- calcNormFactors(tcc, norm.method = "tmm", test.method = "edger",
                       iteration = 3, FDR = 0.1, floorPDEG = 0.05)
plot(tcc, median.lines = TRUE)

