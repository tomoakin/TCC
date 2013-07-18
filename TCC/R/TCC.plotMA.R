
TCC$methods(plotMA = function (FDR = NULL,
                       significance.level = NULL,
                       median.lines = FALSE,
                       floor = 0,
                       groups = NULL,
                       main = NULL,
                       xlab = expression(A == (log[2] * G2 + log[2] * G1 ) / 2),
                       ylab = expression(M == log[2] * G2 - log[2] * G1),
                       xlim = NULL,
                       ylim = NULL,
                       cex = 0.3,
                       pch = 19,
                       col = NULL, 
                       col.tag = NULL, ...) {
    ## set up default arguments.
    gro <- .self$group[, 1]
    gru <- unique(as.vector(gro))
    if (is.null(groups)) {
        groups <- c(gru[1], gru[2])
    }
    if (is.null(col)) {
        if (private$estimated == TRUE) {
            col <- c(1, rep(6, length = length(gru)))
        } else if (private$simulation == TRUE) {
            col <- c(1, 4, 2, 4 + 1:(length(gru)))
        } else {
            col <- rep(1, length = length(gru))
        }
    ##} else {
    ##  if (length(col) != length(gru) + 1) {
    ##    if (length(col) == 1)
    ##      col <- c(col, col)
    ##    col <- c(col[1], rep(col[-1], length = length(gru) - 1))
    ##  }
    }
  
    count.normed <- .self$getNormalizedData()  
    mean.i <- rowMeans(as.matrix(count.normed[, gro == groups[1]]))
    mean.j <- rowMeans(as.matrix(count.normed[, gro == groups[2]]))
    norm.i <- mean(norm.factors[gro == groups[1]])
    norm.j <- mean(norm.factors[gro == groups[2]])
    ma.axes <- .self$.getMACoordinates(mean.i, mean.j, floor)
    filter <- as.logical(mean.i > 0 & mean.j > 0)
    a <- ma.axes$a.value
    m <- ma.axes$m.value
  
    ##m.values <- array(0, dim = c(length(gru), length(gru), nrow(count)))
    ##a.values <- array(0, dim = c(length(gru), length(gru), nrow(count)))
    ## calculate the average of count expression.
    ##mean.exp <- matrix(0, ncol = length(gru), nrow = nrow(count))
    ##mean.norm <- rep(0, length = length(gru))
    ##for (g in 1:length(gru)) {
    ##  mean.exp[, g] <- rowMeans(as.matrix(count.normed[, gro == gru[g]]))
    ##  mean.norm[g] <- mean(norm.factors[gro == gru[g]])
    ##}
    ## calculate m.values and a.values of each combinations of groups.
    ##fig.tils <- length(gru) - 1
    ##if (length(gru) > 2) {
    ##  split.screen(figs = c(fig.tils, fig.tils))
    ##}
    ##global.mar <- par("mar")
    ##global.cex <- par("cex")
    ##for (i in 1:length(gru)) {
    ##  for (j in i:length(gru)) {
    ##    if (i != j) {
    ##      ma.axes <- .self$.getMACoordinates(mean.exp[, i], 
    ##                                         mean.exp[, j], floor)
    ##      filter <- as.logical(mean.exp[, i] > 0 & mean.exp[, j] > 0)
    ##      a.values[i, j, ] <- ma.axes$a.value
    ##      m.values[i, j, ] <- ma.axes$m.value
    ##      a.values[j, i, ] <- ma.axes$a.value
    ##      m.values[j, i, ] <- ma.axes$m.value
    ##      a <- a.values[i, j, ]
    ##      m <- m.values[i, j, ]
    ##      #if (length(gru) > 2) {
    ##      #  screen(fig.tils * (i - 1) + j - 1)
    ##      #}
    ##      #if (length(gru) > 2) {
    ##      #  par(cex = 0.7 * global.cex)
    ##      #}
    ##      #if (length(gru) > 4) {
    ##      #  par(mar = global.mar / (length(gru) - 1))
    ##      #}
    if (is.null(xlim))
        xlim <- c(min(a), max(a))
    if (is.null(ylim))
        ylim <- c(min(m), max(m))
    if (is.null(main))
        main <- "MA plot"
    if (length(gru) > 2) {
        gftitle <- paste(main, " (Group ", groups[1], 
                         " - Group ", groups[2],")", sep = "")
    } else {
        gftitle <- main
    } 
    plot(0, 0, xlim = xlim, ylim = ylim, 
         xlab = xlab, ylab = ylab, type = "n", ...)
    title(main = list(gftitle))
    grid(col = "gray", lty = "dotted")
    ##if (is.null(col.tag)) {
    col.tag.v <- rep(0, length = nrow(count))
    if (private$estimated == FALSE) {
        if (private$simulation == TRUE)
            col.tag.v <- simulation$trueDEG
    } else {
        if ((!is.null(estimatedDEG)) && (length(estimatedDEG != 0))) {
            col.tag.v <- as.numeric(estimatedDEG)
        }
        if (!(is.null(FDR) && is.null(significance.level))) {
            private$stat$q.value <<- stat$q.value
            private$stat$p.value <<- stat$p.value
            col.tag.v <- .self$.exactTest(FDR = FDR, 
                                  significance.level = significance.level)
        }
    }
    ##} else {
    if (is.null(col.tag))
        col.tag <- col.tag.v + 1
    if (length(col.tag) != nrow(count))
        stop("\nTCC::ERROR: The length of col.tag has to be equal to the number of genes.\n")
    for (k in unique(col.tag)) {
        points(a[col.tag == k], m[col.tag == k], 
               col = col[k], pch = pch, cex = cex)
    }
    if (median.lines == TRUE) {
        for (k in unique(col.tag)) {
            if (length(setdiff(gru, groups)) != 0 && k == setdiff(gru, groups))
              next
            med <- median(m[(col.tag == k & filter)])
            lines(c(min(a) + 1, max(a)), c(med, med), col = col[k])
            text(xlim[2], med + 0.5, sprintf("%.3f", med), col = col[k],
                 pos = 2, offset = 0)
        }
    }
    ##    }
    ##  }
    ##}
    ##if (length(gru) > 2) {
    ##  par(mar = global.mar)
    ##  close.screen(all.screens = TRUE)
    ##} 
    ##if (length(gru) == 2) {
    ##  invisible(data.frame(a.value = a.values[i, j, ], 
    ##                       m.value = m.values[i, j, ]))
    invisible(data.frame(a.value = a, m.value = m))
    ##}
})

TCC$methods(.getMACoordinates = function(g1, g2, floor = 0) {
    m <- rep(0, length = nrow(count))
    a <- rep(0, length = nrow(count))
    g1.min.nonZero <- min(g1[g1 > 0])
    g2.min.nonZero <- min(g2[g2 > 0])
    filter <- as.logical(g1 <= floor | g2 <= floor)
    g1[g1 <= floor] <- g1.min.nonZero
    g2[g2 <= floor] <- g2.min.nonZero
    a <- (log2(g1) + log2(g2)) / 2
    m <- log2(g2) - log2(g1)
    a[filter] <- min(a) - 1
    return(list(m.value = m, a.value = a))
})

