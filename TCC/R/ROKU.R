.outval <- function(y, upper.limit) {
    if (all(is.na(y)))
        y <- rep(0, length = length(y))
    z <- y[!is.na(y)]
    
    N <- length(z)
    fN <- floor(N * upper.limit) + 1

    z.sorted <- sort(z)
    z.ordered <- order(z)
    
    df <- matrix(c(rep(1:fN, times = c(fN:1)), j = sequence(fN:1)),
                 ncol = 2, byrow = FALSE)
    n <- N - df[, 2] - df[, 1] + 2
    s <- N - n
    f <- rep(0, length = N)

    if (sd(z) != 0) {
        ssd <- apply(df, 1, function(d, w = z.sorted, N = N) {
           return(sd(w[d[1]:(N + 1 - d[2])]))
        }, z.sorted, N)
        u <- n * log(ssd * sqrt((n - 1) / n)) + 
             sqrt(2) * s * lfactorial(n) / n
        min.u <- min(u)
        d <- t(df[u == min.u, ])[1:2]
        if (d[1] > 1)
            f[z.ordered[1:(d[1] - 1)]] <- -1
        if (d[2] > 1)
            f[z.ordered[(N + 1 - d[2] + 1):N]] <- 1
    } 
    return(replace(y, !is.na(y), f))
}

.entval <- function(y) {
    y <- y[!is.na(y)]
    l <- length(y)
    y <- y[y != 0]
    if (sum(y) <= 0 || sd(y) == 0) {
        return (log2(l))
    } else {
        y.m <- median(y)
        y.u <- (y - y.m) / (5 * median(abs(y - y.m)) + 1e-04)
        y.w <- rep(0, length(y))
        y.i <- abs(y.u) <= 1
        y.w[y.i] <- ((1 - y.u^2)^2)[y.i]
        y.b <- sum(y.w * y) / sum(y.w)
        p <- abs(y - y.b)
        p <- p / sum(p)
        return( - sum(p * log2(p)))
    }
}

ROKU <- function(data, upper.limit = 0.25, sort = FALSE) {
   outliers <- t(apply(t(apply(data, 1, scale)), 1,
                          function (y, upper.limit = upper.limit) {
                          .outval(y, upper.limit = upper.limit)
                     }, upper.limit))
   entropy.score <- apply(data, 1, .entval)
   rank <- rank(entropy.score)
   if (!is.null(colnames(data))) {
       l <- c(colnames(data), "entropy.score", "rank")
   } else {
       l <- c(paste("tissue", 1:ncol(data), sep = "_"), "entropy.score", "rank")
   }
   if (!is.null(rownames(data))) {
       r <- rownames(data)
   } else {
       r <- 1:nrow(data)
   }
   df <- cbind(outliers, entropy.score, rank)
   colnames(df) <- l
   rownames(df) <- r
   if (sort) {
       df <- df[order(rank), ]
   }
   return (df)
}

