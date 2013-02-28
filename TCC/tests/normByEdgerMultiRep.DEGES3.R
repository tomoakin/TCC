if (Sys.getenv("TCC_REAL_TEST") != ""){
  library(TCC)
  tcc <- generateSimulationData(Ngene = 100, PDEG = 0.2, DEG.assign = c(0.9, 0.1))
  tcc <- calcNormFactors(tcc, norm.method = "tmm", test.method = "edger", iteration = 3)
  cat("tcc$norm.factors: ")
  cat(tcc$norm.factors)
  cat("\n")
  cat("tcc$stat$execution.time: ")
  cat(tcc$stat$execution.time)
  cat("\n")
}
