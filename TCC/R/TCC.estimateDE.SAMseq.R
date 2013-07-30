TCC$methods(.testBySamseq = function(...) {

.testBySAMseq.1 = function(samplesize = NULL) {
    c <- round(.self$getNormalizedData())
    s <- samr::SAMseq(x = c, y = .self$group[, 1],
                resp.type = "Two class unpaired",
                nperms = samplesize)
    private$stat$testStat <<- s$samr.obj$tt
    private$stat$p.value <<- rep(NA, length = nrow(.self$count))
    private$stat$q.value <<- rep(NA, length = nrow(.self$count))
    private$stat$rank <<- rank(- abs(s$samr.obj$tt))
}

.testBySAMseq.2 = function(samplesize = NULL) {
    c <- round(.self$getNormalizedData())
    s <- samr::SAMseq(x = c, y = .self$group[, 1],
                resp.type = "Multiclass",
                nperms = samplesize)
    private$stat$testStat <<- s$samr.obj$tt
    private$stat$p.value <<- rep(NA, length = nrow(.self$count))
    private$stat$q.value <<- rep(NA, length = nrow(.self$count))
    private$stat$rank <<- rank(- abs(s$samr.obj$tt))
}

al <- list(...)
if (is.null(al$samplesize)) {
    samplesize <- 10
} else {
    samplesize <- al$samplesize
}
ts <- .self$.testStrategy()
if (ts == 1) {
   .testBySAMseq.1(samplesize = samplesize)
} else if (ts == 2) {
   .testBySAMseq.2(samplesize = samplesize)
} else if (ts == 3) {
   stop()
} else {
   stop()
}
})


