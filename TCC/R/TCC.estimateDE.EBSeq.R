TCC$methods(.testByEbseq = function() {
    g <- .self$group[, 1]
    ug <- unique(g)
    if (length(ug) == 2) {
        suppressMessages(EBout <- EBSeq::EBTest(Data = .self$count,
                     Conditions = as.factor(g),
                     sizeFactors = .self$norm.factors * colSums(.self$count),
                     maxround = 10))
        PP <- EBSeq::GetPPMat(EBout)
        df <- matrix(1, ncol = 2, nrow = nrow(.self$count))
        rownames(df) <- rownames(.self$count)
        df[rownames(PP), 1] <- PP[, 1]
        df[rownames(PP), 2] <- PP[, 2]

        private$stat$p.values <<- df[, 1]
        private$stat$q.values <<- 1 - df[, 2]
        private$stat$rank <<- rank(private$stat$p.values)
    } else {
        gp <- matrix(c(rep(1, length = length(ug)), 1:length(ug)),
                     nrow = 2, byrow = TRUE)
        colnames(gp) <- ug
        rownames(gp) <- c("Pattern1", "Pattern2")
        suppressMessages(MultiOut <- EBSeq::EBMultiTest(.self$count,
                     NgVector = NULL,
                     Conditions = g,
                     AllParti = gp,
                     sizeFactors = .self$norm.factors * colSums(.self$count),
                     maxround = 10))
        PP <- EBSeq::GetMultiPP(MultiOut)
        df <- matrix(1, ncol = 2, nrow = nrow(.self$count))
        rownames(df) <- rownames(.self$count)
        df[rownames(PP$PP), 1] <- PP$PP[, 1]
        df[rownames(PP$PP), 2] <- PP$PP[, 2]

        private$stat$p.values <<- df[, 1]
        private$stat$q.values <<- rep(1, length = nrow(.self$count))
        private$stat$rank <<- rank(private$stat$p.values)
    }
})


