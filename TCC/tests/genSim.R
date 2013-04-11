if (Sys.getenv("TCC_REAL_TEST")!=""){
library(TCC)
  set.seed(1234567)
  tcc <- generateSimulationData(Ngene = 10000, PDEG = 0.3,
                         DEG.assign = c(0.6, 0.2, 0.2),
                         DEG.foldchange = c(3, 10, 6),
                         replicates = c(3, 3, 3))
  cat("dim(tcc$count): ")
  cat(dim(tcc$count))
  cat("\n")
  cat("tcc$group: ")
  print(tcc$group)
  cat("\n")
  cat("tcc$count: ")
  cat("\n")
  print(head(tcc$count))
  png("plotFC-multiGroups-uniform.png", 500, 500)
  plotFCPseudocolor(tcc)
  dev.off()
  png("plotMA-multiGroups-uniform.png", 600, 600)
  plot(tcc)
  dev.off()

}
