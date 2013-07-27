TCC$methods(.testByBayseq = function(samplesize = NULL, cl = NULL,
                                     comparison = NULL) {

.testByBayseq.2 = function(samplesize = NULL, cl = NULL) {
    suppressMessages(d <- new("countData", data = round(.self$count),
             replicates = .self$group[, 1],
             groups = list(NDE = rep(1, length = nrow(.self$group)),
                           DE = .self$group[, 1]),
             libsizes = colSums(.self$count) * .self$norm.factors))
    suppressMessages(suppressMessages(d <- getPriors.NB(d, 
                                       samplesize = samplesize,
                                       estimation = "QL", cl = cl)))
    capture.output(d <- getLikelihoods.NB(d, pET = "BIC", cl = cl))
    stat.bayseq <- topCounts(d, group = "DE", number = nrow(.self$count))
    stat.bayseq <- stat.bayseq[rownames(.self$count), ]
    private$stat$rank <<- rank(- d@posteriors[, "DE"])
    private$stat$likelihood <<- stat.bayseq$Likelihood
    private$stat$p.value <<- 1 - stat.bayseq$Likelihood
    private$stat$q.value <<- stat.bayseq$FDR
    private$estimatedDEG <<- as.numeric(.self$private$stat$rank < 
                                  (nrow(.self$count) * d@estProps[2]))
    private$tbt$estProps <<- d@estProps[2]
}

.testByBayseq.3 = function(samplesize = NULL, cl = NULL,
                                       comparison = NULL) {
    if (is.null(comparison))
        comparison <- colnames(.self$group)[2]
    gs <- .self$group
    gs <- cbind(rep(1, length = nrow(.self$group)), gs)
    colnames(gs)[1] <- "NDE"
    suppressMessages(d <- new("countData", data = round(.self$count),
             replicates = .self$group[, 1],
             groups = gs,
             libsizes = colSums(.self$count) * .self$norm.factors))
    capture.output(suppressMessages(d <- getPriors.NB(d,
                                      samplesize = samplesize,
                                      estimation = "QL", cl = cl)))
    capture.output(suppressMessages(d <- getLikelihoods.NB(d,
                                      pET = "BIC", cl = cl)))
    stat.bayseq <- topCounts(d, group = comparison, number = nrow(.self$count))
    stat.bayseq <- stat.bayseq[rownames(.self$count), ]
    private$stat$rank <<- rank(- d@posteriors[, comparison])
    private$stat$likelihood <<- stat.bayseq$Likelihood
    private$stat$p.value <<- 1 - stat.bayseq$Likelihood
    private$stat$q.value <<- stat.bayseq$FDR
    private$estimatedDEG <<- as.numeric(.self$private$stat$rank < 
                                  (nrow(.self$count) * d@estProps[2]))
    private$tbt$estProps <<- d@estProps[2]
}

ts <- .self$.testStrategy()
if (ts == 1) {
   .testByBayseq.2(samplesize = samplesize, cl = cl)
} else if (ts == 2) {
   .testByBayseq.2(samplesize = samplesize, cl = cl)
} else if (ts == 3) {
   .testByBayseq.3(samplesize = samplesize, cl = cl,
                         comparison = comparison)
} else {
   stop()
}

})


