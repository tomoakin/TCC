
TCC$methods(.testByNbpseq = function() {
    g <- .self$group[, 1]
    ug <- unique(g)
    if (length(ug) == 2) {
        nbp <- nbp.test(counts = .self$count,
                        norm.factors = .self$norm.factors,
                        print.level = 0)
        private$stat$p.values <<- nbp$p.values
        private$stat$q.values <<- nbp$q.values
        private$stat$rank <<- rank(private$stat$p.value)
    }
})



