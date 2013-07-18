
TCC$methods(.wad = function(x, g, AD = FALSE, k = 1) {
    x[x <= 0] <- k
    x <- log2(x)
    ug <- unique(g)
    s <- combn(length(ug), 2, function(ij, x = x, g = g, ug = ug, AD = AD) {
        g1 <- (g == ug[ij[1]])
        g2 <- (g == ug[ij[2]])
        m1 <- rowMeans(as.matrix(x[, g1]))
        m2 <- rowMeans(as.matrix(x[, g2]))
        if (AD == TRUE) {
            x_ave <- abs(m1 - m2) / 2
        } else {
            x_ave <- (m1 + m2) / 2
        }
        weight <- (x_ave - min(x_ave)) / (max(x_ave) - min(x_ave))
        s <- abs((m2 - m1) * weight)
        return (s)
    }, TRUE, x, g, ug, AD)
    s <- apply(s, 1, function(i) {
             return (i[max(abs(i)) == abs(i)])
         })
    return(s)
})

WAD <- function(data, group, k = 1) {
    tcc <- new("TCC", count = data, group = group)
    tcc$.testByWad(k = k)
    return(tcc$private$stat$wad)
}

