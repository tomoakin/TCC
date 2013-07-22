
TCC$methods(.testByWad = function(floor.value) {
    ef <- colSums(count) * norm.factors
    x <- sweep(count, 2, mean(ef) / ef, "*")
    s <- .wad(x = x, group = .self$group[, 1],
              log.scale = TRUE,
              floor.value = floor.value)
    private$stat$rank <<- rank(- abs(s))
    private$stat$wad <<- s
    private$estimatedDEG <<- rep(0, length = nrow(count))
})


