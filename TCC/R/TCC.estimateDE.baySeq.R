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


