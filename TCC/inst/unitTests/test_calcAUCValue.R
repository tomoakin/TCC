test_calcAUCValue <- function() {
    tcc <- simulateReadCounts()
    tcc <- calcNormFactors(tcc)
    tcc <- estimateDE(tcc)
    auc <- calcAUCValue(tcc)

    checkTrue(is.numeric(auc))
}

