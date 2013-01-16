if (Sys.getenv("TCC_REAL_TEST")){
  library(TCC)
  data(hypoData)
  tcc <- generateSimulationData(Ngene = 100, PDEG = 0.2, DEG.assign = c(0.9, 0.1))
#  group <- c(3, 3)
#  tcc <- new("TCC", hypoData, c(3,3))
  tcc$calcNormFactors(norm.method = "tmm", test.method = "bayseq")
  tcc$estimateDE(test.method = "edger", FDR = 0.1)
  result <- getResult(tcc, sort = TRUE)
  print(head(result))
}
