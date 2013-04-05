#  TCC
#  v2.21
#  
#  SUN jianqiang
#  2013-01-03
#
#  Class name:
#    TCC
#
#  Required packagse:
#    edgeR
#    baySeq
#    DESeq
#    ROC

TCC <- setRefClass(
  "TCC",
  fields = list(
    count = "matrix",              # counts data of libraries.
    group = "data.frame",             # group of libraries.
    norm.factors = "numeric",      # normalization factors.
    names = "character",           # gene names
    stat = "list",                 # the result of identify DE genes.
    estimatedDEG = "numeric",      # identified result by identifyDEG().
    simulation = "list",           # the aurgument inputed.
    DEGES = "list",                # detailed informations on DEGES .
    private = "list"         # intermediate data on DEGES process.
  ),

  #  Class Methods.
  methods = list(
    initialize = function(count=NULL, group=NULL, norm.factors=NULL, names=NULL) {
      # If count data setted, fill it to TCC class object.
      if(!is.null(count)){ #guard for copy constructor
        if(is.null(group))
          stop("TCC::ERROR: group must be provided.\n")
        if(!is.data.frame(group))
          .self$group <- data.frame(group = group)
        else
          .self$group <- group
        if (nrow(.self$group) != ncol(count)){
          print(.self$group)
          message(paste("nrow(.self$group):", nrow(.self$group), "ncol(count):", ncol(count)))
          stop("TCC::ERROR: The length of group has to be equal to the columns of count data.\n")
        }
        # Fill count data.
        if(!is.matrix(count)){
          count <<- as.matrix(count)
        }else{
          count <<- count
	}
	# count is a matrix with or without colnames, rownames 
        if (is.null(rownames(count))){
          new_names <- paste("gene_", c(1:nrow(count)), sep="")
          rownames(count) <<- new_names
          names <<- new_names
	} else {
          names <<- rownames(count)
        }
        if (is.null(colnames(count))) {
          g <- as.numeric(table(group))
          colnames(count) <<- paste("G", rep(1:length(g), times = g), "_rep", sequence(g), sep = "")
        } else {
          # if the column is not unique, it occurs error on edgeR.
          colns <- colnames(count)
          if (sum(match(colns, colns)) != sum(1:length(colns))) {
            message("TCC::INFO: Changed the column names of count data to unique.")
            colnames(count) <<- paste(colns, 1:length(colns), sep="_")
          }
        }
        rownames(.self$group) <<- colnames(.self$count)
        # Fill normalization factors.
        if (is.null(norm.factors)) {
          normf <- rep(1, length=ncol(count))
          names(normf) <- colnames(.self$count)
        } else {
          if (length(norm.factors) != ncol(count))
            stop("\nTCC::ERROR: The length of norm.factors has to be equal to the columns of cuont data.\n")
          normf <- norm.factors
        }
        norm.factors <<- normf
        private <<- list(estimated = FALSE, simulation = FALSE, normalized = FALSE)
      }
    },
    #/**
    # * THE METHODS OF CALCULATE NORMALIZATION FACTORS.
    # */
    #  TMM normalization. (edgeR)
    .normByTmm = function (count) {
      #if (!("edgeR" %in% loadedNamespaces()))
      #  library(edgeR)
      suppressMessages(d <- edgeR::DGEList(counts=round(count), group=group[[1]]))
      suppressMessages(d <- edgeR::calcNormFactors(d))
      normf <- d$samples$norm.factors
      names(normf) <- colnames(.self$count)
      normf
    },
    #  DESeq normalization. (DESeq) // Each cols is divided by the genomic means of the rows.
    .normByDeseq = function(count) {
      suppressMessages(d <- newCountDataSet(countData=count, conditions=group[[1]]))
      suppressMessages(d <- estimateSizeFactors(d))
      return(sizeFactors(d) / colSums(count))
    }
  )
)
    #  Non-parametric exact test by baySeq.
TCC$methods(    .testByBayseq = function(samplesize, processors, comparison = NULL, cl = NULL) {
      if (is.null(comparison))
          comparison <- colnames(group)[1]
      d <- new("countData",
          data = round(count),
          replicates = group[[1]],
          groups = c(list(NDE = rep(1, length = nrow(group))),
          as.list(group)),
          libsizes = colSums(count) * norm.factors)
      suppressMessages(d <- getPriors.NB(d, samplesize=samplesize, estimation="QL", cl=cl))
      capture.output(suppressMessages(d <- getLikelihoods.NB(d, pET="BIC", cl=cl)))
      stat.bayseq <- topCounts(d, group=comparison, number=nrow(count))
      stat.bayseq <- stat.bayseq[rownames(count), ]
      private$stat$rank <<- rank(- d@posteriors[, comparison])
      # calculate p.value and q.value from likelihood values.
      private$stat$likelihood <<- stat.bayseq$Likelihood
      private$stat$p.value <<- 1 - stat.bayseq$Likelihood
      private$stat$q.value <<- stat.bayseq$FDR
      private$estimatedDEG <<- as.numeric(private$stat$rank < (nrow(count) * d@estProps[2]))
      private$tbt$estProps <<- d@estProps[2]
    }
)

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
      mean.exp <- matrix(0, ncol=length(table(group)), nrow=nrow(count))
      replicates = table(group)
      for (g in 1:length(replicates)) {
        mean.exp[, g] <- log2(rowMeans(as.matrix(count.normed[, replicates == g])))
      }
      ngroup <- length(replicates)
      for (i in 1:(ngroup-1)) {
        for (j in (i+1):ngroup) {
            log2ration <- mean.exp[, j] - mean.exp[, i]
            deg.flg[(deg.flg > 0) & (log2ration < 0)] <- i
            deg.flg[(deg.flg > 0) & (log2ration > 0)] <- j
        }
      }
      return (deg.flg)
    }
)
 
    # exact test.
TCC$methods(estimateDE = function (test.method=NULL,
                           FDR = NULL,
                           significance.level = NULL,
                           samplesize=10000,
                           comparison=NULL,
                           cl=NULL) {
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
        "bayseq" = .self$.testByBayseq(samplesize, comparison=comparison, cl=cl),
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
    }
)
TCC$methods(getNormalizedCount = function () {
      effective.libsizes <- colSums(count) * norm.factors
      return (sweep(count, 2, mean(effective.libsizes) / effective.libsizes, "*"))
    }
)
TCC$methods(.getMACoordinates = function(g1, g2, floor = 0) {
       m <- rep(0, length = nrow(count))
       a <- rep(0, length = nrow(count))
       g1.min.nonZero <- min(g1[g1 > 0])
       g2.min.nonZero <- min(g2[g2 > 0])
       filter <- as.logical(g1 <= floor | g2 <= floor)
       g1[g1 <= floor] <- g1.min.nonZero
       g2[g2 <= floor] <- g2.min.nonZero
       print(min(g1))
       a <- (log2(g1) + log2(g2)) / 2
       m <- log2(g2) - log2(g1)
       a[filter] <- min(a) - 1
       return(list(m.value = m, a.value = a))
    }
)

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
      ngroup <- length(table(group))
      if (is.null(col)) {
        if (private$estimated == TRUE) {
          col <- c(1, rep(6, length= ngroup))
        } else if (private$simulation == TRUE) {
          col <- c(1, 4, 2, 4 + 1:(ngroup - 3))
        } else {
          col <- rep(1, length=ngroup)
        }
      } else {
        if (length(col) != ngroup + 1) {
          if (length(col) == 1)
            col <- c(col, col)
          col <- c(col[1], rep(col[-1], length=ngroup - 1))
        }
      }
      m.values <- array(0, dim=c(ngroup, ngroup, nrow(count)))
      a.values <- array(0, dim=c(ngroup, ngroup, nrow(count)))

      # calculate the average of count expression.
      count.normed <- getNormalizedCount()
      mean.exp <- matrix(0, ncol=ngroup, nrow=nrow(count))
      mean.norm <- rep(0, length=ngroup)
      replicates <- table(group)
      for (g in 1:ngroup) {
        mean.exp[, g] <- rowMeans(as.matrix(count.normed[, group == g]))
        mean.norm[g] <- mean(norm.factors[group == g])
      }
      # calculate m.values and a.values of each combinations of groups.
      fig.tils <- ngroup - 1
      if (ngroup > 2) {
        split.screen(figs = c(fig.tils, fig.tils))
      }
      global.mar <- par("mar")
      global.cex <- par("cex")
      for (i in 1:(ngroup - 1)) {
        for (j in (i+1):ngroup) {
            ma.axes <- .self$.getMACoordinates(mean.exp[, i], mean.exp[, j], floor)
            filter <- as.logical(mean.exp[, i] > 0 & mean.exp[, j] > 0)
            a.values[i, j, ] <- ma.axes$a.value
            m.values[i, j, ] <- ma.axes$m.value
            a.values[j, i, ] <- ma.axes$a.value
            m.values[j, i, ] <- ma.axes$m.value
            a <- a.values[i, j, ]
            m <- m.values[i, j, ]
            if (ngroup > 2) {
              screen(fig.tils * (i - 1) + j - 1)
            }
            if (ngroup > 2) {
              par(cex = 0.7 * global.cex)
            }
            if (ngroup > 4) {
              par(mar = global.mar / (ngroup - 1))
            }
            if (is.null(xlim))
              xlim <- c(min(a), max(a))
            if (is.null(ylim))
              ylim <- c(min(m), max(m))
            if (is.null(main))
              main <- "MA plot"
            if (ngroup > 2) {
               gftitle <- paste(main, " (Group ", i, " - Group ", j,")", sep="")
            } else {
               gftitle <- main
            } 
            plot(0, 0, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, type="n", ...)
            title(main = list(gftitle))
            grid(col="gray", lty="dotted")

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
      if (ngroup > 2) {
        par(mar = global.mar)
        close.screen(all.screens = TRUE)
      } 
      if (ngroup == 2) {
        invisible(data.frame(m.value = m.values[1, 2, ], a.value = a.values[1, 2, ]))
      }
    }
)
