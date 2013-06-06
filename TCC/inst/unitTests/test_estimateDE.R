test_estimateDE <- function() {
    data(hypoData)
    group <- c(1, 1, 1, 2, 2, 2)
    tcc <- new("TCC", hypoData, group)
    tcc <- calcNormFactors(tcc)
    tcc <- estimateDE(tcc, test.method = "edger")
    tcc <- estimateDE(tcc, test.method = "deseq")
    tcc <- estimateDE(tcc, test.method = "bayseq", samplesize = 10)

    data(hypoData_mg)
    group <- c(1, 1, 1, 2, 2, 2, 3, 3, 3)
    tcc <- new("TCC", hypoData_mg, group)
    tcc <- calcNormFactors(tcc)
    desgin <- model.matrix(~ 0 * factor(group))
    tcc <- estimateDE(tcc, test.method = "edger", design = design)
    fit1 <- count ~ condition
    fit0 <- count ~ 1
    tcc <- estimateDE(tcc, test.method = "deseq", fit1 = fit1, fit0 = fit0)
    tcc <- estimateDE(tcc, test.method = "bayseq", samplesize = 10)
}

