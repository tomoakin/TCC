if (Sys.getenv("TCC_REAL_TEST")!=""){
  library(TCC)
  tcc <- generateSimulationData(Ngene = 100, PDEG = 0.2, DEG.assign = c(0.9, 0.1))
  tcc <- calcNormFactors(tcc)
  tcc <- estimateDE(tcc, test.method = "bayseq", FDR = 0.1, samplesize = 100)
  result <- getResult(tcc, sort = TRUE)
  print(head(result))
  table(tcc$estimatedDEG) 
  png("plotMA-testByBayseq.png", 600, 500)
  plot(tcc)
  dev.off()
}
