TCC <- setRefClass(
  "TCC",
  fields = list(
    count = "matrix",              # counts data of libraries.
    gene_id = "character",         # gene names
    group = "data.frame",          # groups, libraries, conditions.
    norm.factors = "numeric",      # normalization factors.
    stat = "list",                 # the result of identify DE genes.
    estimatedDEG = "numeric",      # identified result by identifyDEG().
    simulation = "list",           # the aurgument inputed.
    DEGES = "list",                # detailed informations about DEGES .
    private = "list",              # intermediate data on DEGES process.
    debug = "list"
  ),

  #  Class Methods.
  methods = list(
    initialize = function(count = NULL, group = NULL,
                          norm.factors = NULL, gene_id = NULL) {
      if (!is.null(count)) {
        #if (is.null(group) && is.null(replicates))
        #  stop("TCC::ERROR: The group or replicates must be provided.\n")
        # Set up group data.
        if (is.null(group)) {
          stop("TCC::ERROR: The group or replicates must be provided.\n")
          #.self$group <<- data.frame(group = rep(1:length(replicates), times = replicates))
        } else {
          if (!is.data.frame(group)) 
            .self$group <<- data.frame(group = group)
          else 
            .self$group <<- group
        }
        # Set up count data.
        if (!is.matrix(count))
          .self$count <<- as.matrix(count)
        else
          .self$count <<- count
        # Set up names.
        if (is.null(rownames(.self$count))) {
          .self$gene_id <<- paste("gene_", c(1:nrow(count)), sep = "")
          rownames(.self$count) <<- paste("gene_", c(1:nrow(count)), sep = "")
        } else {
          .self$gene_id <<- rownames(count)
        }
        if (is.null(colnames(.self$count))) {
          g <- as.numeric(table(group))
          colnames(.self$count) <<- paste("G", rep(1:length(g), times = g), "_rep", sequence(g), sep = "")
        } else {
          colnm <- colnames(count)
          if (sum(match(colnm, colnm)) != sum(1:length(colnm))) {
            message("TCC::INFO: Changed the column names of count data to unique.")
            colnames(count) <<- paste(colnm, 1:length(colnm), sep = "_")
          }
        }
        rownames(.self$group) <<- colnames(.self$count)
        # Set up normlization factors if it was given.
        if (is.null(norm.factors)) {
          .self$norm.factors <<- rep(1, length = ncol(count))
        } else {
          if (length(norm.factors) != ncol(count))
            stop("\nTCC::ERROR: The length of norm.factors has to be equal to the columns of cuont data.\n")
          .self$norm.factors <<- norm.factors
        }
        names(norm.factors) <<- norm.factors
        # Set private argument.
        private$estimated <<- FALSE
        private$simulation <<- FALSE
        private$normalized <<- FALSE
      }
    }
  )
)

#/**
# * THE METHODS OF CALCULATE NORMALIZATION FACTORS.
# */
TCC$methods(.normByTmm = function(count){
  suppressMessages(d <- edgeR::DGEList(counts = round(count), group = group[, 1]))
  suppressMessages(d <- edgeR::calcNormFactors(d))
  normf <- d$samples$norm.factors
  names(normf) <- colnames(.self$count)
  return(normf)
})
TCC$methods(.normByDeseq = function(count){
  if (ncol(group) == 1)
    suppressMessages(d <- newCountDataSet(countData = round(count), conditions = group[, 1]))
  else 
    suppressMessages(d <- newCountDataSet(countData = round(count), conditions = group))
  suppressMessages(d <- estimateSizeFactors(d))
  return(sizeFactors(d) / colSums(count))
})



#/**
# * THE METHODS OF INDENTIFY DE GENES.
# */
#  Parametric exact test by edgeR.
TCC$methods(.testByEdger = function(design = NULL, coef = NULL, contrast = NULL, dispersion = NULL){
  if (is.null(contrast) && (is.null(coef))) {
    suppressMessages(d <- edgeR::DGEList(counts = round(count), group = group[, 1]))
    suppressMessages(d <- edgeR::calcNormFactors(d))
    d$samples$norm.factors <- norm.factors
    if (min(table(group[, 1])) > 1) {
      suppressMessages(d <- edgeR::estimateCommonDisp(d))
      suppressMessages(d <- edgeR::estimateTagwiseDisp(d))
    }
    if (is.null(dispersion)) {
      suppressMessages(d <- edgeR::exactTest(d))
    } else {
      suppressMessages(d <- edgeR::exactTest(d, dispersion = dispersion))
    }
    if (!is.null(d$table$PValue)) {
      private$stat$p.value <<- d$table$PValue
    } else {
      private$stat$p.value <<- d$table$p.value
    }
    private$stat$rank <<- rank(private$stat$p.value)
    private$stat$q.value <<- p.adjust(private$stat$p.value, method = "BH")
  } else {
    if (is.null(design))
      stop("TCC::ERROR: Need the design matrix for GLM.")
    # GLM test.
    suppressMessages(d <- edgeR::DGEList(counts = round(count), group = group[, 1]))
    suppressMessages(d <- edgeR::calcNormFactors(d))
    d$samples$norm.factors <- norm.factors
    suppressMessages(d <- edgeR::estimateGLMCommonDisp(d, design))
    suppressMessages(d <- edgeR::estimateGLMTrendedDisp(d, design))
    suppressMessages(d <- edgeR::estimateGLMTagwiseDisp(d, design))
    suppressMessages(fit <- edgeR::glmFit(d, design))
    suppressMessages(lrt <- edgeR::glmLRT(fit, coef = coef, contrast = contrast))
    s <- topTags(lrt, n = nrow(count))
    s <- s$table[rownames(count), ]
    private$stat$p.value <<- s$PValue
    private$stat$rank <<- rank(private$stat$p.value)
    private$stat$q.value <<- s$FDR
  }
})
#  Parametric exact test by DESeq.
TCC$methods(.testByDeseq = function(fit1 = NULL, fit0 = NULL, comparison = NULL){
  if (is.null(comparison))
    comparison <- colnames(group)[1]
  if (ncol(group) == 1)
    suppressMessages(d <- newCountDataSet(countData = round(count), conditions = group[, 1]))
  else 
    suppressMessages(d <- newCountDataSet(countData = round(count), conditions = group))
  sizeFactors(d) <- norm.factors * colSums(count)
  if (ncol(group) == 1 && min(as.numeric(table(group[, 1]))) == 1) { # single group and single replicate
    e <- try(suppressMessages(d <- estimateDispersions(d, method = "blind", sharingMode = "fit-only")), silent = TRUE)
    if (class(e) == "try-error") {
      message("TCC::WARN: 'estimateDispersions' with sharingMode=\"fit-only\" in DESeq could not be performed.")
      message("TCC::WARN: 'estimateDispersions' with sharingMode=\"local\" in DESeq was used instead.")
      suppressMessages(d <- estimateDispersions(d, fitType = "local"))
    }
  } else { # otherwise conditions
    # try default
    e <- try(suppressMessages(d <- estimateDispersions(d)), silent = TRUE)
    # try blind method
    if (class(e) == "try-error") {
      message("TCC::WARN: 'estimateDispersions' with method=\"pooled\" in DESeq could not be performed.")
      message("TCC::WARN: 'estimateDispersions' with method=\"blind\" in DESeq was used instead.")
      suppressMessages(e <- (d <- estimateDispersions(d, method = "blind", sharingMode = "fit-only")))
      # try local mode
      if (class(e) == "try-error") {
        message("TCC::WARN: 'estimateDispersions' with sharingMode=\"fit-only\" in DESeq could not be performed.")
        message("TCC::WARN: 'estimateDispersions' with sharingMode=\"local\" in DESeq was used instead.")
        suppressMessages(d <- estimateDispersions(d, fitType = "local"))
      }
    }
  }
  # classic or GLM
  if (is.null(fit1) && is.null(fit0)) {
    unique.group <- sort(unique(group[, comparison]))
    suppressMessages(d <- nbinomTest(d, unique.group[1], unique.group[2]))
    d$pval[is.na(d$pval)] <- 1
    d$padj[is.na(d$padj)] <- 1
    private$stat$p.value <<- d$pval
    private$stat$q.value <<- d$padj
    private$stat$rank <<- rank(d$pval)
  } else {
    if (is.null(fit0))
      stop("TCC::ERROR: Need the formula('fit0') to create reduced model regresses for GLM.")
    if (is.null(fit1))
      stop("TCC::ERROR: Need the formula('fit1') to create full model regresses for GLM.")
    capture.output(fit0 <- fitNbinomGLMs(d, fit0))
    capture.output(fit1 <- fitNbinomGLMs(d, fit1))
    private$stat$p.value <<- nbinomGLMTest(fit1, fit0)
    private$stat$p.value[is.na(private$stat$p.value)] <<- 1
    private$stat$q.value <<- p.adjust(private$stat$p.value, method = "BH")
    private$stat$rank <<- rank(private$stat$p.value)
  }
})
#  Non-parametric exact test by baySeq.
TCC$methods(.testByBayseq = function(samplesize = NULL, cl = NULL, comparison = NULL){
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
  suppressMessages(d <- getPriors.NB(d, samplesize = samplesize, estimation = "QL", cl = cl))
  capture.output(suppressMessages(d <- getLikelihoods.NB(d, pET = "BIC", cl = cl)))
  stat.bayseq <- topCounts(d, group = comparison, number = nrow(count))
  stat.bayseq <- stat.bayseq[rownames(count), ]
  private$stat$rank <<- rank(- d@posteriors[, comparison])
  # calculate p.value and q.value from likelihood values.
  private$stat$likelihood <<- stat.bayseq$Likelihood
  private$stat$p.value <<- 1 - stat.bayseq$Likelihood
  private$stat$q.value <<- stat.bayseq$FDR
  private$estimatedDEG <<- as.numeric(private$stat$rank < (nrow(count) * d@estProps[2]))
  private$tbt$estProps <<- d@estProps[2]
})

#  calculate normalization factors.
TCC$methods(calcNormFactors = function(norm.method = NULL,
                                test.method = NULL,
                                iteration = 1,
                                FDR = NULL,
                                floorPDEG = NULL,
                                dispersion = NULL,
                                design = NULL,
                                contrast = NULL, coef = NULL,
                                fit0 = NULL, fit1 = NULL,
                                comparison = NULL,
                                samplesize = 10000,
                                cl = NULL, increment = FALSE){
  if ((increment == FALSE) || 
      (increment == TRUE && private$normalized == FALSE)) {
    DEGES$iteration <<- 0
    private$debug <<- list()
    private$debug$potentialDEG <<- list()
    private$debug$potentialPval <<- list()
  }
  ex.time <- proc.time()
  if (is.null(norm.method)) {
    if ((ncol(group) == 1) && (min(as.numeric(table(group))) == 1)) 
      norm.method = "deseq"
    else 
      norm.method = "edger"
  }
  if (is.null(test.method)) {
    if ((ncol(group) == 1) && (min(as.numeric(table(group))) == 1)) 
      test.method = "deseq"
    else 
      test.method = "edger"
  }
  if (norm.method == "edger")
    norm.method <- "tmm" 
  if (test.method != "bayseq" && is.null(FDR))
    FDR <- 0.1
  if (test.method != "bayseq" && is.null(floorPDEG)) 
    floorPDEG <- 0.05
  if (iteration) {
    if (is.logical(iteration))
      iteration <- 1
    if (iteration < 0 && 100 < iteration)
      stop("TCC::ERROR: The iteration must be given within the limits of from '0' to '100'.")
    message(paste("TCC::INFO: Calculating normalization factors using DEGES"))
    message(paste("TCC::INFO: (iDEGES pipeline :", norm.method, 
                  "- [", test.method, "-", norm.method, "] X", iteration + DEGES$iteration, ")"))
    DEGES$protocol <<- paste(norm.method, "- [", test.method, "-", norm.method, "] X", iteration + DEGES$iteration)
  } else {
    message(paste("TCC::INFO: Calculating normalization factors using", norm.method, "..."))
    DEGES$protocol <<- norm.method
  }
  # DEGES strategy STEP 1. (First normalization)
  if ((increment == FALSE) || 
      (increment == TRUE && private$normalized == FALSE)) {
    norm.factors <<- switch(norm.method,
      "tmm" = .self$.normByTmm(count),
      "deseq" = .self$.normByDeseq(count),
      stop(paste("\nTCC::ERROR: The normalization method of ", norm.method, " doesn't supported.\n"))
    )
  }
  norm.factors <<- norm.factors / mean(norm.factors)
  DEGES$threshold <<- data.frame(type = "Unused", input = 0, PDEG = 0)
  #  if DEGES not NULL then start DEGES strategy.
  if (iteration) {
    #  if iteration > 0 then change to iterate DEGES strategy.
    for (i in 1:iteration) {
      # DEGES strategy  STEP 2. (exact test and remove differential expression genes.)
      private$stat <<- list()
      switch(test.method,
        "edger" = .self$.testByEdger(design = design, coef = coef, contrast = contrast, dispersion = dispersion),
        "deseq" = .self$.testByDeseq(fit1 = fit1, fit0 = fit0, comparison = comparison),
        "bayseq" = .self$.testByBayseq(samplesize = samplesize, cl = cl, comparison = comparison),
        stop(paste("\nTCC::ERROR: The identifying method of ", test.method, " doesn't supported.\n"))
      )
      # Remove the DEG from original count data.
      deg.flg.FDR <- .self$.exactTest(FDR = FDR)
      deg.flg.floorPDEG <- as.numeric(rank(private$stat$p.value, ties.method = "min") <= nrow(count) * floorPDEG)
      if (is.null(floorPDEG) && is.null(FDR)) {
        # use TbT
        deg.flg <- deg.flg.FDR
        DEGES$threshold$type <<- "TbT"
        DEGES$threshold$input <<- private$tbt$estProps
        DEGES$threshold$PDEG <<- sum(deg.flg) / length(deg.flg)
      } else {
        # use FDR
        deg.flg <- deg.flg.FDR
        DEGES$threshold$type <<- "FDR"
        DEGES$threshold$input <<- FDR
        DEGES$threshold$PDEG <<- sum(deg.flg) / length(deg.flg)
        if ((!is.null(floorPDEG)) && (sum(deg.flg != 0) < sum(deg.flg.floorPDEG != 0))) {
        # use floorPDEG
          deg.flg <- deg.flg.floorPDEG
          DEGES$threshold$type <<- "floorPDEG"
          DEGES$threshold$input <<- floorPDEG
          DEGES$threshold$PDEG <<- sum(deg.flg) / length(deg.flg)
        }
      }
      count.ndeg <- count[(deg.flg == 0), ]
      if (nrow(count.ndeg) == 0) {
        message ("TCC::INFO: No non-DE genes after eliminate DE genes. stop DEGES strategy.")
        break
      }
      private$debug$potentialDEG[[DEGES$iteration + 1]] <<- deg.flg
      private$debug$potentialThreshold[[DEGES$iteration + 1]] <<- DEGES$threshold
      private$debug$potentialPval[[DEGES$iteration + 1]] <<- private$stat$p.value
      # DEGES strategy STEP 3. (Second normalization)
      norm.factors <<- switch(norm.method,
        "tmm" = .self$.normByTmm(count.ndeg),
        "deseq" = .self$.normByDeseq(count.ndeg)
      )
      norm.factors <<- norm.factors * colSums(count.ndeg) / colSums(count)
      norm.factors <<- norm.factors / mean(norm.factors)
      DEGES$iteration <<- DEGES$iteration + 1
    }
    DEGES$potentialDEG <<- deg.flg
  }
  message("TCC::INFO: Done.")
  DEGES$execution.time <<- proc.time() - ex.time
  private$normalized <<- TRUE
})

TCC$methods(.exactTest = function (FDR = NULL, significance.level = NULL) {
  deg.flg <- rep(0, length = nrow(count))
  if (!is.null(significance.level)) {
    deg.flg <- as.numeric(private$stat$p.value < significance.level)
  } else if (!is.null(FDR)) {
    deg.flg <- as.numeric(private$stat$q.value < FDR)
  } else {
    deg.flg <- private$estimatedDEG #TbT
  }
  return (deg.flg)
})
 
    # exact test.
TCC$methods(estimateDE = function (test.method = NULL,
                           FDR = NULL,
                           significance.level = NULL,
                           dispersion = NULL,
                           fit0 = NULL, fit1 = NULL,
                           design = NULL,
                           contrast = NULL, coef = NULL,
                           comparison = NULL,
                           samplesize = 10000,
                           cl = NULL) {
  if (is.null(test.method)) {
    if ((ncol(group) == 1) && (min(as.numeric(table(group))) == 1)) 
      test.method = "deseq"
    else 
      test.method = "edger"
  }
  if (test.method != "bayseq" && is.null(FDR) && is.null(significance.level)) 
    FDR <- 0.1
  message(paste("TCC::INFO: Identifying DE genes using", test.method, "..."))
  # calculate statistics values related DE gene.
  private$stat <<- list()
  switch(test.method,
    "edger" = .self$.testByEdger(design = design, coef = coef, contrast = contrast, dispersion = dispersion),
    "deseq" = .self$.testByDeseq(fit1 = fit1, fit0 = fit0, comparison = comparison),
    "bayseq" = .self$.testByBayseq(samplesize = samplesize, cl = cl, comparison = comparison),
    stop(paste("\nTCC::ERROR: The identifying method of ", test.method, " doesn't supported.\n"))
  )
  # identify DE genes with the results of exact test.
  estimatedDEG <<- .self$.exactTest(FDR = FDR, significance.level = significance.level)
  if (!is.null(private$stat$likelihood))
    stat$likelihood <<- private$stat$likelihood
  if (!is.null(private$stat$p.value))
    stat$p.value <<- private$stat$p.value
  if (!is.null(private$stat$q.value))
    stat$q.value <<- private$stat$q.value
  if (!is.null(private$stat$rank))
    stat$rank <<- private$stat$rank
  private$estimated <<- TRUE
  message("TCC::INFO: Done.")
})
TCC$methods(getNormalizedCount = function () {
  effective.libsizes <- colSums(count) * norm.factors
  return (sweep(count, 2, mean(effective.libsizes) / effective.libsizes, "*"))
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

# plot M-A plotting.
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
  # set up default arguments.
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
  } else {
    if (length(col) != length(gru) + 1) {
      if (length(col) == 1)
        col <- c(col, col)
      col <- c(col[1], rep(col[-1], length = length(gru) - 1))
    }
  }
  
  count.normed <- getNormalizedCount()  
  mean.i <- rowMeans(as.matrix(count.normed[, gro == groups[1]]))
  mean.j <- rowMeans(as.matrix(count.normed[, gro == groups[2]]))
  norm.i <- mean(norm.factors[gro == groups[1]])
  norm.j <- mean(norm.factors[gro == groups[2]])
  ma.axes <- .self$.getMACoordinates(mean.i, mean.j, floor)
  filter <- as.logical(mean.i > 0 & mean.j > 0)
  a <- ma.axes$a.value
  m <- ma.axes$m.value
  
  #m.values <- array(0, dim = c(length(gru), length(gru), nrow(count)))
  #a.values <- array(0, dim = c(length(gru), length(gru), nrow(count)))
  # calculate the average of count expression.
  #mean.exp <- matrix(0, ncol = length(gru), nrow = nrow(count))
  #mean.norm <- rep(0, length = length(gru))
  #for (g in 1:length(gru)) {
  #  mean.exp[, g] <- rowMeans(as.matrix(count.normed[, gro == gru[g]]))
  #  mean.norm[g] <- mean(norm.factors[gro == gru[g]])
  #}
  # calculate m.values and a.values of each combinations of groups.
  #fig.tils <- length(gru) - 1
  #if (length(gru) > 2) {
  #  split.screen(figs = c(fig.tils, fig.tils))
  #}
  #global.mar <- par("mar")
  #global.cex <- par("cex")
  #for (i in 1:length(gru)) {
  #  for (j in i:length(gru)) {
  #    if (i != j) {
  #      ma.axes <- .self$.getMACoordinates(mean.exp[, i], mean.exp[, j], floor)
  #      filter <- as.logical(mean.exp[, i] > 0 & mean.exp[, j] > 0)
  #      a.values[i, j, ] <- ma.axes$a.value
  #      m.values[i, j, ] <- ma.axes$m.value
  #      a.values[j, i, ] <- ma.axes$a.value
  #      m.values[j, i, ] <- ma.axes$m.value
  #      a <- a.values[i, j, ]
  #      m <- m.values[i, j, ]
        #if (length(gru) > 2) {
        #  screen(fig.tils * (i - 1) + j - 1)
        #}
        #if (length(gru) > 2) {
        #  par(cex = 0.7 * global.cex)
        #}
        #if (length(gru) > 4) {
        #  par(mar = global.mar / (length(gru) - 1))
        #}
        if (is.null(xlim))
          xlim <- c(min(a), max(a))
        if (is.null(ylim))
          ylim <- c(min(m), max(m))
        if (is.null(main))
          main <- "MA plot"
        if (length(gru) > 2) {
           gftitle <- paste(main, " (Group ", groups[1], " - Group ", groups[2],")", sep = "")
        } else {
           gftitle <- main
        } 
        plot(0, 0, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, type = "n", ...)
        title(main = list(gftitle))
        grid(col = "gray", lty = "dotted")
        #if (is.null(col.tag)) {
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
              col.tag.v <- .self$.exactTest(FDR = FDR, significance.level = significance.level)
            }
          }
		#} else {
        if (is.null(col.tag))
          col.tag <- col.tag.v
        if (length(col.tag) != nrow(count))
            stop("\nTCC::ERROR: The length of col.tag has to be equal to the number of genes.\n")
        for (k in 0:max(col.tag)) {
          points(a[col.tag == k], m[col.tag == k], col = col[k + 1], pch = pch, cex = cex)
        }
        if (median.lines == TRUE) {
          for (k in unique(col.tag.v)) {
            if (length(setdiff(gru, groups)) != 0 && k == setdiff(gru, groups))
              next
            med <- median(m[(col.tag == k & filter)])
            lines(c(min(a) + 1, max(a)), c(med, med), col = col[k + 1])
            text(xlim[2], med + 0.5, sprintf("%.3f", med), col = col[k + 1],
                 pos = 2, offset = 0)
          }
        }
  #    }
  #  }
  #}
  #if (length(gru) > 2) {
  #  par(mar = global.mar)
  #  close.screen(all.screens = TRUE)
  #} 
  #if (length(gru) == 2) {
  #  invisible(data.frame(a.value = a.values[i, j, ], m.value = m.values[i, j, ]))
    invisible(data.frame(a.value = a, m.value = m))
  #}
})
