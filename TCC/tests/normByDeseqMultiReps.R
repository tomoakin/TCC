if (Sys.getenv("TCC_REAL_TEST")!=""){
  library(TCC)
  tcc <- simulateReadCounts(Ngene = 100, PDEG = 0.2, DEG.assign = c(0.9, 0.1))
  tcc <- calcNormFactors(tcc, norm.method = "deseq", test.method = "deseq",
                          FDR = 0.1, floorPDEG = 0.05)
  cat("tcc$norm.factors: ")
  cat(tcc$norm.factors)
  cat("\n")
  cat("tcc$DEGES$execution.time: ")
  cat(tcc$DEGES$execution.time)
  cat("\n")
}
