TCC$methods(.testByNoiseq = function() {
    g <- .self$group[, 1]
    ug <- unique(g)
    if (length(ug) == 2) {
        x <- .self$getNormalizedData()
        gl <- data.frame(group = g)
        nd <- NOISeq::readData(x, gl)
        suppressMessages(nr <- NOISeq::noiseq(nd, k = 0.5,
                         norm = "n", replicates = "biological",
                         factor = "group", conditions = ug))
        private$stat$p.values <<- nr@results[[1]]$prob
        private$stat$q.values <<- rep(1, length = nrow(.self$count))
        private$stat$rank <<- rank(- private$stat$p.values)
    }
})


