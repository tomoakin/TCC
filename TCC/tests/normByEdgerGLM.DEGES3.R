if (Sys.getenv("TCC_REAL_TEST") != ""){
  library(TCC)
  group <- c(1, 1, 1, 2, 2, 2, 3, 3, 3)
  replicates <- c(3, 3, 3)
  design <- model.matrix(~ factor(group))
  tcc <- simulateReadCounts(Ngene = 100, PDEG = 0.2, DEG.assign = c(0.6, 0.2, 0.2),
                                replicates = replicates)
  tcc <- calcNormFactors(tcc, norm.method = "tmm", test.method = "edger", 
                         design = design, coef = 2:3, iteration = 3)
  cat("tcc$norm.factors: ")
  cat(tcc$norm.factors)
  cat("\n")

  design <- model.matrix(~ 0 + factor(group))
  tcc <- calcNormFactors(tcc, norm.method = "tmm", test.method = "edger", 
                         design = design, contrast = c(-1, 0, 1), iteration = 3)
  cat("tcc$norm.factors: ")
  cat(tcc$norm.factors)
  cat("\n")
}
