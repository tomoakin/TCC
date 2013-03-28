if (Sys.getenv("TCC_REAL_TEST") != ""){
  library(TCC)
  replicates <- c(3, 3, 3)
  fit1 <- count ~ condition
  fit0 <- count ~ 1
  tcc <- generateSimulationData(Ngene = 100, PDEG = 0.2, DEG.assign = c(0.6, 0.2, 0.2),
                                replicates = replicates)
  tcc <- calcNormFactors(tcc, norm.method = "deseq", test.method = "deseq", 
                         fit1 = fit1, fit0 = fit0, iteration = 3)
  cat("tcc$norm.factors: ")
  cat(tcc$norm.factors)
  cat("\n")
}
