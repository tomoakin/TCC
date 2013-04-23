library(TCC)
set.seed(137427)
tcc <- simulateReadCounts(Ngene = 100, PDEG = 0.3,
                       DEG.assign = c(0.6, 0.2, 0.2),
                       DEG.foldchange = c(3, 10, 6),
                       replicates = c(2, 4, 3))
cat("dim(tcc$count): ")
cat(dim(tcc$count))
cat("\n")
cat("tcc$group$group: ")
cat(tcc$group$group)
cat("\n")
cat("tcc$count:\n")
print(head(tcc$count))
plotFCPseudocolor(tcc)
plot(tcc)

tcc <- simulateReadCounts(Ngene = 100, PDEG = 0.3,
                       DEG.assign = c(0.6, 0.2, 0.2),
                       DEG.foldchange = c(3, 10, 6),
                       replicates = c(2, 4, 3))
cat("dim(tcc$count): ")
cat(dim(tcc$count))
cat("\n")
cat("tcc$group$group: ")
cat(tcc$group$group)
cat("\n")
cat("tcc$count:\n")
print(head(tcc$count))
plotFCPseudocolor(tcc)
plot(tcc)

tcc <- simulateReadCounts(Ngene = 200, PDEG = 0.30,
                              DEG.assign = c(0.85, 0.15),
                              DEG.foldchange = list(min = c(1.2, 1.2),
                                                    shape = c(2.0, 2.0),
                                                    scale = c(0.5, 0.5)),
                              replicates = c(2, 2))
tcc <- calcNormFactors(tcc, norm.method = "tmm", test.method = "edger",
                       iteration = 3, FDR = 0.1, floorPDEG = 0.05)
plot(tcc, median.lines = TRUE)

