TCC <- setRefClass(
  "TCC",
  fields = list(
    count = "matrix",              # counts data of libraries.
    names = "character",           # gene names
    group = "numeric",             # group of libraries.
    replicates = "numeric",        # group of libraries.
    norm.factors = "numeric",      # normalization factors.
    stat = "list",                 # the result of identify DE genes.
    estimatedDEG = "numeric",      # identified result by identifyDEG().
    simulation = "list",           # the aurgument inputed.
    private = "list"         # intermediate data on DEGES process.
  ),

  #  Class Methods.
  methods = list(
    initialize = function(count = NULL, group = NULL, replicates = NULL, norm.factors = NULL, names = NULL) {
      # If count data setted, fill it to TCC class object.
      if (!is.null(count)) {
        if(!is.null(group)){
          if (sum(group) != ncol(count))
            stop("TCC::ERROR: The sum of group has to be equal to the columns of count data.\n")
          replicates <<- rep(1:length(group), times = group)
          group <<- group
        }else if(!is.null(replicates)){
          if (length(replicates) != ncol(count))
            stop("TCC::ERROR: The length of replicates has to be equal to the columns of count data.\n")
          group <<- rep(0, length = max(replicates))
          for (i in 1:max(replicates)) {
            group[i] <<- sum(replicates == i)
          }
          replicates <<- replicates
        }else{
          stop("TCC::ERROR: group or replicate must be provided.\n")
        }
        # Fill count data.
        if(!is.matrix(count)){
          count <<- as.matrix(count)
        }else{
          count <<- count
	}
	# count is a matrix with or without colnames, rownames 
        if (is.null(rownames(count))){
          names <<- paste("gene_", c(1:nrow(count)), sep="")
          rownames(count) <<- names
          names <<- names
	} else {
          names <<- rownames(count)
        }
        if (is.null(colnames(count))) {
          colnames(count) <<- paste("G", rep(1:length(group), times = group), "_rep", sequence(group), sep = "")
        } else {
          # if the column is not unique, it occurs error on edgeR.
          colns <- colnames(count)
          if (sum(match(colns, colns)) != sum(1:length(colns))) {
            message("TCC::INFO: Changed the column names of count data to unique.")
            colnames(count) <<- paste(colns, 1:length(colns), sep = "_")
          }
        }
        # Fill normalization factors.
        if (is.null(norm.factors)) {
          normf <- rep(1, length = ncol(count))
	  names(normf) <- colnames(count)
        } else {
          if (length(norm.factors) != ncol(count))
            stop("\nTCC::ERROR: The length of norm.factors has to be equal to the columns of cuont data.\n")
          normf <- norm.factors
        }
        norm.factors <<- normf
      }
      private$estimated <<- FALSE
      private$simulation <<- FALSE
    },
    #/**
    # * THE METHODS OF CALCULATE NORMALIZATION FACTORS.
    # */
    #  TMM normalization. (edgeR)
    .normByTmm = function (count) {
      #if (!("edgeR" %in% loadedNamespaces()))
      #  library(edgeR)
      suppressMessages(d <- edgeR::DGEList(counts = count, group = replicates))
      suppressMessages(d <- edgeR::calcNormFactors(d))
      normf <- d$samples$norm.factors
      names(normf) <- colnames(.self$count)
      normf
    },
    #  DESeq normalization. (DESeq) // Each cols is divided by the genomic means of the rows.
    .normByDeseq = function(count) {
      suppressMessages(d <- newCountDataSet(countData = count, conditions = replicates))
      suppressMessages(d <- estimateSizeFactors(d))
      return(sizeFactors(d) / colSums(count))
    }
  )
)
#/**
# * THE METHODS OF INDENTIFY DE GENES.
# */
#  Parametric exact test by edgeR.
TCC$methods(.testByEdger = function(){
  suppressMessages(d <- edgeR::DGEList(counts = count, group = replicates))
  suppressMessages(d <- edgeR::calcNormFactors(d))
  d$samples$norm.factors <- norm.factors
  suppressMessages(d <- edgeR::estimateCommonDisp(d))
  suppressMessages(d <- edgeR::estimateTagwiseDisp(d))
  suppressMessages(d <- edgeR::exactTest(d))
  if (!is.null(d$table$PValue)) {
    private$stat$p.value <<- d$table$PValue
  } else {
    private$stat$p.value <<- d$table$p.value
  }
  private$stat$rank <<- rank(private$stat$p.value)
  private$stat$q.value <<- p.adjust(private$stat$p.value, method = "BH")
})
#  Parametric exact test by DESeq.
TCC$methods(.testByDeseq = function(){
  suppressMessages(d <- newCountDataSet(countData = count, conditions = replicates))
  sizeFactors(d) <- norm.factors * colSums(count)
  if (ncol(count) > 2) {
    e <- try(suppressMessages(d <- estimateDispersions(d)), silent = TRUE)
    if (class(e) == "try-error") {
      message("TCC::WARN: An Error occurs when execute 'estimateDispersions' in DESeq.")
      message("TCC::WARN: Changed 'fitType' to 'local' of 'estiamteDispersions'.")
      suppressMessages(d <- estimateDispersions(d, fitType = "local"))
    }
  } else {
    e <- try(suppressMessages(d <- estimateDispersions(d, method = "blind", sharingMode = "fit-only")), silent = TRUE)
    if (class(e) == "try-error") {
      message("TCC::WARN: An Error occurs when execute 'estimateDispersions' in DESeq.")
      message("TCC::WARN: Changed 'fitType' to 'local' of 'estiamteDispersions'.")
      suppressMessages(d <- estimateDispersions(d, method = "blind", sharingMode = "fit-only", fitType = "local"))
    }
  }
  suppressMessages(d <- nbinomTest(d, 1, 2))
  d$pval[is.na(d$pval)] <- 1
  d$padj[is.na(d$padj)] <- 1
  private$stat$p.value <<- d$pval
  private$stat$q.value <<- d$padj
  private$stat$rank <<- rank(d$pval)
})
#  Non-parametric exact test by baySeq.
TCC$methods(.testByBayseq = function(samplesize, processors){
  cl <- NULL
  if (!is.null(processors)) {
    if(!("snow" %in% loadedNamespaces()))
      library(snow)
    if (is.numeric(processors)) {
      cl <- makeSOCKcluster(rep("localhost", length = processors))
    } else {
      cl <- processors
    }
  }
  suppressMessages(d <- new("countData", data = as.matrix(count), 
      replicates = replicates, 
      groups = list(NDE = rep(1, length = length(replicates)), DE = replicates), 
      libsizes = colSums(count) * norm.factors))
  suppressMessages(d <- getPriors.NB(d, samplesize = samplesize, estimation = "QL", cl = cl))
  capture.output(suppressMessages(d <- getLikelihoods.NB(d, pET = "BIC", cl = cl)))
  stat.bayseq <- topCounts(d, group = "DE", number = nrow(count))
  stat.bayseq <- stat.bayseq[rownames(count), ]
  private$stat$rank <<- rank(- d@posteriors[, "DE"])
  # calculate p.value and q.value from likelihood values.
  private$stat$likelihood <<- stat.bayseq$Likelihood
  private$stat$p.value <<- 1 - stat.bayseq$Likelihood
  private$stat$q.value <<- stat.bayseq$FDR
  private$estimatedDEG <<- as.numeric(private$stat$rank < (nrow(count) * d@estProps[2]))
})

#  calculate normalization factors.
TCC$methods(calcNormFactors = function(norm.method = NULL,
                                test.method = NULL,
                                iteration = 1,
                                FDR = NULL,
                                floorPDEG = NULL,
                                samplesize = 10000,
                                processors = NULL){
  if (is.null(norm.method)) {
    if (min(.self$group) == 1) {
      norm.method = "deseq"
    } else {
      norm.method = "edger"
    }
  }
  if (is.null(test.method)) {
    if (ncol(count) < 4) {
      test.method = "deseq"
    } else {
      test.method = "edger"
    }
  }
  if (norm.method == "edger")
    norm.method <- "tmm" 
  if (test.method != "bayseq" && is.null(FDR)) {
    FDR <- 0.1
  }
  if (test.method != "bayseq" && is.null(floorPDEG)) {
    floorPDEG <- 0.05
  }
  if (iteration) {
    if (is.logical(iteration))
      iteration <- 1
    message(paste("TCC::INFO: Calculating normalization factors using DEGES"))
    message(paste("TCC::INFO: (iDEGES pipeline :", norm.method, "- [", test.method, "-", norm.method, "] X", iteration, ")"))
  } else {
    message(paste("TCC::INFO: Calculating normalization factors using", norm.method, "..."))
  }
  # DEGES strategy STEP 1. (First normalization)
  norm.factors <<- switch(norm.method,
    "tmm" = .self$.normByTmm(count),
    "deseq" = .self$.normByDeseq(count),
    stop(paste("\nTCC::ERROR: The normalization method of ", norm.method, " doesn't supported.\n"))
  )
  private$DEGES.threshold.type <<- 0
  #  if DEGES not NULL then start DEGES strategy.
  if (iteration) {
    #  if iteration > 0 then change to iterate DEGES strategy.
    for (i in 1:iteration) {
      # DEGES strategy  STEP 2. (exact test and remove differential expression genes.)
      private$stat <<- list()
      switch(test.method,
        "edger" = .self$.testByEdger(),
        "deseq" = .self$.testByDeseq(),
        "bayseq" = .self$.testByBayseq(samplesize, processors),
        stop(paste("\nTCC::ERROR: The identifying method of ", test.method, " doesn't supported.\n"))
      )
      # Remove the DEG from original count data.
      deg.flg.FDR <- .self$.exactTest(FDR = FDR)
      deg.flg.floorPDEG <- as.numeric(rank(private$stat$p.value, ties.method = "min") <= nrow(count) * floorPDEG)
      if (is.null(FDR)) {
        deg.flg <- deg.flg.FDR
        private$DEGES.threshold.type <<- 5 + sum(private$estimatedDEG != 0) / length(private$estimatedDEG)
      } else {
        deg.flg <- deg.flg.FDR
        private$DEGES.threshold.type <<- 3 + FDR
      }
      # super threshold.
      if ((!is.null(floorPDEG)) && (sum(deg.flg != 0) < sum(deg.flg.floorPDEG != 0))) {
        deg.flg <- deg.flg.floorPDEG
        private$DEGES.threshold.type <<- 1 + floorPDEG
      }
      count.ndeg <- count[(deg.flg == 0), ]
      if (nrow(count.ndeg) == 0) {
        message ("TCC::INFO: No non-DE genes after eliminate DE genes. stop DEGES strategy.")
        break
      }
      # DEGES strategy STEP 3. (Second normalization)
      norm.factors <<- switch(norm.method,
        "tmm" = .self$.normByTmm(count.ndeg),
        "deseq" = .self$.normByDeseq(count.ndeg)
      )
      norm.factors <<- norm.factors * colSums(count.ndeg) / colSums(count)
      norm.factors <<- norm.factors / mean(norm.factors)
    }
    private$DEGES.potentialDEG <<- deg.flg
  } else {
    norm.factors <<- norm.factors / mean(norm.factors)
  }
  message("TCC::INFO: Done.")
})

TCC$methods(.exactTest = function (FDR = NULL, significance.level = NULL) {
  deg.flg <- rep(0, length = nrow(count))
  if (!is.null(significance.level)) {
    deg.flg <- as.numeric(private$stat$p.value < significance.level)
    private$estimate.type <<- 1 + significance.level
  } else if (!is.null(FDR)) {
    deg.flg <- as.numeric(private$stat$q.value < FDR)
    private$estimate.type <<- 3 + FDR
  } else {
    # Only for TbT
    deg.flg <- private$estimatedDEG
    private$estimate.type <<- 5 + sum(private$estimatedDEG != 0) / length(private$estimatedDEG)
  }
  # decide group of DEG.
  count.normed <- .self$getNormalizedCount()
  mean.exp <- matrix(0, ncol = length(group), nrow = nrow(count))
  for (g in 1:length(group)) {
    mean.exp[, g] <- log2(rowMeans(as.matrix(count.normed[, replicates == g])))
  }
  for (i in 1:length(group)) {
    for (j in i:length(group)) {
      if (i != j) {
        log2ration <- mean.exp[, j] - mean.exp[, i]
        deg.flg[(deg.flg > 0) & (log2ration < 0)] <- i
        deg.flg[(deg.flg > 0) & (log2ration > 0)] <- j
      }
    }
  }
  return (deg.flg)
})
 
    # exact test.
TCC$methods(estimateDE = function (test.method = NULL,
                           FDR = NULL,
                           significance.level = NULL,
                           samplesize = 10000,
                           processors = NULL) {
  if (is.null(test.method)) {
    if (ncol(count) < 4) {
       test.method = "deseq"
    } else {
       test.method = "edger"
    }
  }
  if (test.method != "bayseq" && is.null(FDR) && is.null(significance.level)) {
    FDR <- 0.1
  }
  message(paste("TCC::INFO: Identifying DE genes using", test.method, "..."))
  # calculate statistics values related DE gene.
  private$stat <<- list()
  switch(test.method,
    "edger" = .self$.testByEdger(),
    "deseq" = .self$.testByDeseq(),
    "bayseq" = .self$.testByBayseq(samplesize, processors),
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
                       main = NULL,
                       xlab = expression(A == (log[2] * G2 + log[2] * G1 ) / 2),
                       ylab = expression(M == log[2] * G2 - log[2] * G1),
                       xlim = NULL,
                       ylim = NULL,
                       cex = 0.3,
                       pch = 19,
                       col = NULL, ...) {
  # set up default arguments.
  if (is.null(col)) {
    if (private$estimated == TRUE) {
      col <- c(1, rep(6, length = length(group)))
    } else if (private$simulation == TRUE) {
      col <- c(1, 4, 2, 4 + 1:(length(group) - 3))
    } else {
      col <- rep(1, length = length(group))
    }
  } else {
    if (length(col) != length(group) + 1) {
      if (length(col) == 1)
        col <- c(col, col)
      col <- c(col[1], rep(col[-1], length = length(group) - 1))
    }
  }
  m.values <- array(0, dim = c(length(group), length(group), nrow(count)))
  a.values <- array(0, dim = c(length(group), length(group), nrow(count)))

  # calculate the average of count expression.
  count.normed <- getNormalizedCount()
  mean.exp <- matrix(0, ncol = length(group), nrow = nrow(count))
  mean.norm <- rep(0, length = length(group))
  for (g in 1:length(group)) {
    mean.exp[, g] <- rowMeans(as.matrix(count.normed[, replicates == g]))
    mean.norm[g] <- mean(norm.factors[replicates == g])
  }
  # calculate m.values and a.values of each combinations of groups.
  fig.tils <- length(group) - 1
  if (length(group) > 2) {
    split.screen(figs = c(fig.tils, fig.tils))
  }
  global.mar <- par("mar")
  global.cex <- par("cex")
  for (i in 1:length(group)) {
    for (j in i:length(group)) {
      if (i != j) {
        ma.axes <- .self$.getMACoordinates(mean.exp[, i], mean.exp[, j], floor)
        filter <- as.logical(mean.exp[, i] > 0 & mean.exp[, j] > 0)
        a.values[i, j, ] <- ma.axes$a.value
        m.values[i, j, ] <- ma.axes$m.value
        a.values[j, i, ] <- ma.axes$a.value
        m.values[j, i, ] <- ma.axes$m.value
        a <- a.values[i, j, ]
        m <- m.values[i, j, ]
        if (length(group) > 2) {
          screen(fig.tils * (i - 1) + j - 1)
        }
        if (length(group) > 2) {
          par(cex = 0.7 * global.cex)
        }
        if (length(group) > 4) {
          par(mar = global.mar / (length(group) - 1))
        }
        if (is.null(xlim))
          xlim <- c(min(a), max(a))
        if (is.null(ylim))
          ylim <- c(min(m), max(m))
        if (is.null(main))
          main <- "MA plot"
        if (length(group) > 2) {
           gftitle <- paste(main, " (Group ", i, " - Group ", j,")", sep = "")
        } else {
           gftitle <- main
        } 
        plot(0, 0, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, type = "n", ...)
        title(main = list(gftitle))
        grid(col = "gray", lty = "dotted")

        col.tag <- rep(0, length = nrow(count))
        if (private$estimated == FALSE) {
          if (private$simulation == TRUE)
            col.tag <- simulation$trueDEG
        } else {
          if ((!is.null(estimatedDEG)) && (length(estimatedDEG != 0))) {
            col.tag <- as.numeric(estimatedDEG)
          } 
          if (!(is.null(FDR) && is.null(significance.level))) {
            private$stat$q.value <<- stat$q.value
            private$stat$p.value <<- stat$p.value
            col.tag <- .self$.exactTest(FDR = FDR, significance.level = significance.level)
          }
        }
        for (k in 0:max(col.tag)) {
          points(a[col.tag == k], m[col.tag == k], col = col[k + 1], pch = pch, cex = cex)
        }
        if (median.lines == TRUE) {
          for (k in c(0, i, j)) {
            med <- median(m.values[i, j, (col.tag == k & filter)])
            lines(c(min(a) + 1, max(a)), c(med, med), col = col[k + 1])
            text(xlim[2], med + 0.5, sprintf("%.3f", med), col = col[k + 1],
                 pos = 2, offset = 0)
          }
        }
      }
    }
  }
  if (length(group) > 2) {
    par(mar = global.mar)
    close.screen(all.screens = TRUE)
  } 
  if (length(group) == 2) {
    invisible(data.frame(m.value = m.values[1, 2, ], a.value = a.values[1, 2, ]))
  }
})
