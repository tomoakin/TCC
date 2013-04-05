#/**
# * THE METHODS OF INDENTIFY DE GENES.
# */
#  Parametric exact test by edgeR.
TCC$methods(.testByEdger = function(design = NULL, coef = NULL, contrast = NULL, dispersion = NULL){
  if (is.null(contrast) && (is.null(coef))) {
    d <- edgeR::DGEList(counts = round(count), group = group[, 1])
    d <- edgeR::calcNormFactors(d)
    d$samples$norm.factors <- norm.factors
    if (min(table(group[, 1])) > 1) {
      d <- edgeR::estimateCommonDisp(d)
      d <- edgeR::estimateTagwiseDisp(d)
    }
    if (is.null(dispersion)) {
      (d <- edgeR::exactTest(d))
    } else {
      (d <- edgeR::exactTest(d, dispersion = dispersion))
    }
    if (!is.null(d$table$PValue)) {
      private$stat$p.value <<- d$table$PValue
    } else {
      private$stat$p.value <<- d$table$p.value
    }
    private$stat$rank <<- rank(private$stat$p.value)
    private$stat$q.value <<- p.adjust(private$stat$p.value, method = "BH")
  } else {
    if (is.null(design))
      stop("TCC::ERROR: Need the design matrix for GLM.")
    # GLM test.
    suppressMessages(d <- edgeR::DGEList(counts = round(count), group = group[, 1]))
    suppressMessages(d <- edgeR::calcNormFactors(d))
    d$samples$norm.factors <- norm.factors
    suppressMessages(d <- edgeR::estimateGLMCommonDisp(d, design))
    suppressMessages(d <- edgeR::estimateGLMTrendedDisp(d, design))
    suppressMessages(d <- edgeR::estimateGLMTagwiseDisp(d, design))
    suppressMessages(fit <- edgeR::glmFit(d, design))
    suppressMessages(lrt <- edgeR::glmLRT(fit, coef = coef, contrast = contrast))
    s <- topTags(lrt, n = nrow(count))
    s <- s$table[rownames(count), ]
    private$stat$p.value <<- s$PValue
    private$stat$rank <<- rank(private$stat$p.value)
    private$stat$q.value <<- s$FDR
  }
})
