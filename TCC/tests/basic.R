if (Sys.getenv("TCC_REAL_TEST")!=""){
  library(TCC)
  data(hypoData)
  tcc <- simulateReadCounts(Ngene = 100, PDEG = 0.2, DEG.assign = c(0.9, 0.1))
  print(tcc)
  tcc$calcNormFactors(norm.method = "tmm", test.method = "bayseq")
  print(tcc)
  tcc$estimateDE(test.method = "edger", FDR = 0.1)
  print(tcc)
  result <- getResult(tcc, sort = TRUE)
  print(head(result))

  print(names(tcc))
  print(head(tcc[c("gene_1", "gene_2", "gene_3")]))
  print(length(tcc))
  print(head(tcc[c(1,2,3)]))
}
