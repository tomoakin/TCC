TCC$methods(.testByEdger = function(design = NULL, coef = NULL, 
                                    contrast = NULL, dispersion = NULL){
    if (is.null(contrast) && (is.null(coef))) {
        suppressMessages(d <- edgeR::DGEList(counts = round(count), 
                                             group = group[, 1]))
        suppressMessages(d <- edgeR::calcNormFactors(d))
        d$samples$norm.factors <- norm.factors
        if (min(table(group[, 1])) > 1) {
            suppressMessages(d <- edgeR::estimateCommonDisp(d))
            suppressMessages(d <- edgeR::estimateTagwiseDisp(d))
        }
        if (is.null(dispersion)) {
            suppressMessages(d <- edgeR::exactTest(d))
        } else {
            suppressMessages(d <- edgeR::exactTest(d, dispersion = dispersion))
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
        ## GLM test.
        suppressMessages(d <- edgeR::DGEList(counts = round(count), 
                                             group = group[, 1]))
        suppressMessages(d <- edgeR::calcNormFactors(d))
        d$samples$norm.factors <- norm.factors
        suppressMessages(d <- edgeR::estimateGLMCommonDisp(d, design))
        suppressMessages(d <- edgeR::estimateGLMTrendedDisp(d, design))
        suppressMessages(d <- edgeR::estimateGLMTagwiseDisp(d, design))
        suppressMessages(fit <- edgeR::glmFit(d, design))
        suppressMessages(lrt <- edgeR::glmLRT(fit, coef = coef, 
                                              contrast = contrast))
        s <- topTags(lrt, n = nrow(count))
        s <- s$table[rownames(count), ]
        private$stat$p.value <<- s$PValue
        private$stat$rank <<- rank(private$stat$p.value)
        private$stat$q.value <<- s$FDR
    }
})

TCC$methods(.testByDeseq = function(fit1 = NULL, fit0 = NULL, 
                                    comparison = NULL){
    if (is.null(comparison))
        comparison <- colnames(group)[1]
    if (ncol(group) == 1) {
        suppressMessages(d <- newCountDataSet(countData = round(count), 
                                              conditions = group[, 1]))
    } else {
        suppressMessages(d <- newCountDataSet(countData = round(count), 
                                              conditions = group))
    }
    sizeFactors(d) <- norm.factors * colSums(count)
    if (ncol(group) == 1 && min(as.numeric(table(group[, 1]))) == 1) {
        ## single group and single replicate
        e <- try(suppressMessages(d <- estimateDispersions(d, 
                                           method = "blind", 
                                           sharingMode = "fit-only")), 
                 silent = TRUE)
        if (class(e) == "try-error") {
             message("TCC::WARN: 'estimateDispersions' with sharingMode=\"fit-only\" in DESeq could not be performed.")
             message("TCC::WARN: 'estimateDispersions' with sharingMode=\"local\" in DESeq was used instead.")
             suppressMessages(d <- estimateDispersions(d, 
                                       method = "blind", 
                                       sharingMode = "fit-only", 
                                       fitType = "local"))
        }
    } else { 
        ## otherwise conditions
        ## try default
        e <- try(suppressMessages(d <- estimateDispersions(d)), silent = TRUE)
        ## try blind method
        if (class(e) == "try-error") {
            message("TCC::WARN: 'estimateDispersions' with method=\"pooled\" in DESeq could not be performed.")
            message("TCC::WARN: 'estimateDispersions' with method=\"blind\" in DESeq was used instead.")
            e <- try(suppressMessages(d <- estimateDispersions(d, 
                                               method = "blind", 
                                               sharingMode = "fit-only")), 
                     silent = TRUE)
            ## try local mode
            if (class(e) == "try-error") {
                message("TCC::WARN: 'estimateDispersions' with sharingMode=\"fit-only\" in DESeq could not be performed.")
                message("TCC::WARN: 'estimateDispersions' with sharingMode=\"local\" in DESeq was used instead.")
                suppressMessages(d <- estimateDispersions(d, 
                                          method = "blind", 
                                          sharingMode = "fit-only", 
                                          fitType = "local"))
            }
        }
    }
    ## classic or GLM
    if (is.null(fit1) && is.null(fit0)) {
        unique.group <- sort(unique(group[, comparison]))
        suppressMessages(d <- nbinomTest(d, unique.group[1], unique.group[2]))
        d$pval[is.na(d$pval)] <- 1
        d$padj[is.na(d$padj)] <- 1
        private$stat$p.value <<- d$pval
        private$stat$q.value <<- d$padj
        private$stat$rank <<- rank(d$pval)
    } else {
        if (is.null(fit0))
            stop("TCC::ERROR: Need the formula('fit0') to create reduced model regresses for GLM.")
        if (is.null(fit1))
            stop("TCC::ERROR: Need the formula('fit1') to create full model regresses for GLM.")
        capture.output(fit0 <- fitNbinomGLMs(d, fit0))
        capture.output(fit1 <- fitNbinomGLMs(d, fit1))
        private$stat$p.value <<- nbinomGLMTest(fit1, fit0)
        private$stat$p.value[is.na(private$stat$p.value)] <<- 1
        private$stat$q.value <<- p.adjust(private$stat$p.value, method = "BH")
        private$stat$rank <<- rank(private$stat$p.value)
    }
})

TCC$methods(.testByBayseq = function(samplesize = NULL, 
                                     cl = NULL, comparison = NULL){
    if (is.null(comparison)) 
        comparison <- colnames(group)[1]
    suppressMessages(d <- new("countData",
                              data = round(count), 
                              replicates = group[, 1], 
                              groups = c(
                                list(NDE = rep(1, length = nrow(group))),
                                as.list(group)
                              ),
                              libsizes = colSums(count) * norm.factors))
    suppressMessages(d <- getPriors.NB(d, samplesize = samplesize, 
                                       estimation = "QL", cl = cl))
    capture.output(suppressMessages(d <- getLikelihoods.NB(d, pET = "BIC", 
                                                           cl = cl)))
    stat.bayseq <- topCounts(d, group = comparison, number = nrow(count))
    stat.bayseq <- stat.bayseq[rownames(count), ]
    private$stat$rank <<- rank(- d@posteriors[, comparison])
    ## calculate p.value and q.value from likelihood values.
    private$stat$likelihood <<- stat.bayseq$Likelihood
    private$stat$p.value <<- 1 - stat.bayseq$Likelihood
    private$stat$q.value <<- stat.bayseq$FDR
    private$estimatedDEG <<- as.numeric(private$stat$rank < 
                                        (nrow(count) * d@estProps[2]))
    private$tbt$estProps <<- d@estProps[2]
})

TCC$methods(.testByWad = function(k) {
    ef <- colSums(count) * norm.factors
    x <- sweep(count, 2, mean(ef) / ef, "*")
    s <- .self$.wad(x = x, g = .self$group[, 1], k = k)
    private$stat$rank <<- rank(- abs(s))
    private$stat$wad <<- s
    private$estimatedDEG <<- rep(0, length = nrow(count))
})

TCC$methods(.exactTest = function (FDR = NULL, significance.level = NULL,
                                   PDEG = NULL) {
    deg.flg <- rep(0, length = nrow(count))
    if (!is.null(significance.level)) {
        deg.flg <- as.numeric(private$stat$p.value < significance.level)
    } else if (!is.null(FDR)) {
        deg.flg <- as.numeric(private$stat$q.value < FDR)
    } else if (!is.null(PDEG)) {
        deg.flg <- as.numeric(private$stat$rank <= nrow(count) * PDEG)
    } else {
        deg.flg <- private$estimatedDEG #TbT
    }
    return (deg.flg)
})
 

TCC$methods(estimateDE = function (test.method = NULL,
                                   FDR = NULL, PDEG = NULL,
                                   significance.level = NULL,
                                   dispersion = NULL,
                                   fit0 = NULL, fit1 = NULL,
                                   design = NULL,
                                   contrast = NULL, coef = NULL,
                                   comparison = NULL,
                                   samplesize = 10000,
                                   k = 1,
                                   cl = NULL) {
    if (is.null(test.method)) {
        if ((ncol(group) == 1) && (min(as.numeric(table(group))) == 1)) 
            test.method = "deseq"
        else 
            test.method = "edger"
    }
    if (test.method == "wad") {
        PDEG <- 0.1
    } else if (test.method != "bayseq" && is.null(FDR) && 
               is.null(significance.level)) {
        FDR <- 0.1
    }
    message(paste("TCC::INFO: Identifying DE genes using", test.method, "..."))
    ## calculate statistics values related DE gene.
    private$stat <<- list()
    stat <<- list()
    switch(test.method,
           "edger" = .self$.testByEdger(design = design, 
                                        coef = coef, 
                                        contrast = contrast,  
                                        dispersion = dispersion),
           "deseq" = .self$.testByDeseq(fit1 = fit1, 
                                        fit0 = fit0, 
                                        comparison = comparison),
           "bayseq" = .self$.testByBayseq(samplesize = samplesize, 
                                          cl = cl, 
                                          comparison = comparison),
           "wad" = .self$.testByWad(k = k),
           stop(paste("\nTCC::ERROR: The identifying method of ", 
                      test.method, " doesn't supported.\n"))
    )
    ## identify DE genes with the results of exact test.
    estimatedDEG <<- .self$.exactTest(FDR = FDR, 
                                      significance.level = significance.level,
                                      PDEG = PDEG)
    if (!is.null(private$stat$likelihood))
        stat$likelihood <<- private$stat$likelihood
    if (!is.null(private$stat$p.value))
        stat$p.value <<- private$stat$p.value
    if (!is.null(private$stat$q.value))
        stat$q.value <<- private$stat$q.value
    if (!is.null(private$stat$wad)) 
        stat$wad <<- private$stat$wad
    if (!is.null(private$stat$rank)) 
        stat$rank <<- private$stat$rank
    private$estimated <<- TRUE
    message("TCC::INFO: Done.")
})

