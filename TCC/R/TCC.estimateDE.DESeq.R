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


