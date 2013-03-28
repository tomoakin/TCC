#TCC <- function(count = NULL, conditions = NULL, norm.factors = NULL, names = NULL) {
#  tcc <- new("TCC", count = count, conditions = conditions,
#             norm.factors = norm.factors, names = names)
#  return (tcc)
#}


# calcNormFactors
# calculate normalization factors with TCC class tcc.
setGeneric(name = "calcNormFactors", def = function(tcc, ...) tcc)
setMethod(
  f = "calcNormFactors",
  signature(tcc = "DGEList"),
  definition = function(tcc, ...) {
    return(edgeR::calcNormFactors(tcc, ...))
  }
)
setMethod(
  f = "calcNormFactors",
  signature(tcc = "TCC"),
  definition = function(tcc, norm.method=NULL, test.method=NULL, iteration=TRUE,
                        FDR=NULL, floorPDEG=NULL, dispersion=NULL,
                        design=NULL,contrast=NULL, coef=NULL,
                        fit0=NULL, fit1=NULL,  comparison = NULL,
                        samplesize=10000, cl=NULL) {
      obj <- tcc$copy()
      obj$calcNormFactors(norm.method=norm.method, test.method=test.method, iteration=iteration,
                      FDR=FDR, floorPDEG=floorPDEG, dispersion=dispersion,
                      fit0=fit0, fit1=fit1, design=design, contrast=contrast, coef=coef,
                      comparison=comparison,samplesize=samplesize, cl=cl)
      return(obj)
    }
)

# estimateDE
# the method is for estimating DEGs.
estimateDE <- function(tcc, test.method=NULL, FDR=NULL, dispersion=NULL,
                       fit0=NULL, fit1=NULL, design = NULL, contrast=NULL, coef=NULL,
                       comparison=NULL,samplesize=10000, cl=NULL) {
  obj <- tcc$copy()
  obj$estimateDE(test.method=test.method, FDR=FDR, dispersion=dispersion,
                 fit0=fit0, fit1=fit1, design=design, contrast=contrast, coef=coef,
                 comparison=comparison, samplesize=samplesize, cl=cl)
  return(obj)
}

# plot
# plot MA-plot with TCC class tcc.
plot.TCC <- function(x, FDR=NULL, median.lines = FALSE, floor=0, groups=NULL, main=NULL, 
                    xlab = expression(A == (log[2] * G2 + log[2] * G1 ) / 2),
                    ylab = expression(M == log[2] * G2 - log[2] * G1),
                    xlim = NULL, ylim = NULL, cex = 0.3, pch = 19, col = NULL, col.tag = NULL,...) {
      invisible(x$plotMA(FDR=FDR, median.lines=median.lines, floor=floor, groups=groups, main=main, xlab=xlab, ylab=ylab,
               xlim=xlim, ylim=ylim, cex=cex, pch=pch, col=col, col.tag=col.tag,...))
}

# getResult
# get p-value, FDR or the axes of MA-plot as data.frame.
getResult <- function(tcc, sort = FALSE, floor = 0) {
  if (length(tcc$stat) == 0)
    stop("\nTCC::ERROR: There are no statistics in stat fields of TCC class tcc. Execute TCC.estiamteDE for calculating them.\n")
  gru <- unique(tcc$group[, 1])
  if ((length(gru) == 2) && (ncol(tcc$group) == 1)) {
    count.normed <- tcc$getNormalizedCount()
    mean.exp <- matrix(0, ncol=length(gru), nrow=nrow(tcc$count))
    for (g in 1:length(gru))
      mean.exp[, g] <- rowMeans(as.matrix(count.normed[, tcc$group[, 1] == g]))
    ma.axes <- tcc$.getMACoordinates(mean.exp[, 1], mean.exp[, 2], floor)
    result.df <- data.frame(
      id = rownames(tcc$count),
      a.value = ma.axes$a.value, m.value = ma.axes$m.value,
      p.value = tcc$stat$p.value, q.value = tcc$stat$q.value,
      rank = tcc$stat$rank, estimatedDEG = tcc$estimatedDEG
    )
  } else {
    result.df <- data.frame(
      id = rownames(tcc$count),
	  a.value = rep(NA, length = nrow(tcc$count)), m.value = rep(NA, length = nrow(tcc$count)),
      p.value = tcc$stat$p.value, q.value = tcc$stat$q.value,
      rank = tcc$stat$rank, estimatedDEG = tcc$estimatedDEG
    )
  }
  if (sort)
    result.df <- result.df[order(result.df$rank), ]
  return (result.df)
}

# filterData
# remove the low count data.
filterLowCountGenes <- function(tcc, low.count = 0) {
  obj <- tcc$copy()
  gru <- unique(obj$group[, 1])
  filters <- matrix(0, ncol=length(gru), nrow=nrow(obj$count)) 
  for (i in 1:length(gru)) {
    filters[, i] <- as.numeric(rowSums(
                      as.matrix(obj$count[, (obj$group[, 1] == gru[i])])
                    ) <= low.count)
  }
  left.tag <- as.logical(rowSums(filters) != length(gru))
  obj$count <- obj$count[left.tag, ]
  if (!is.null(obj$simulation$trueDEG) && length(obj$simulation$trueDEG) != 0)
    obj$simulation$trueDEG <- obj$simulation$trueDEG[left.tag]
  if (!is.null(obj$estimatedDEG) && length(obj$estimatedDEG) != 0)
    obj$estimatedDEG <- obj$estimatedDEG[left.tag]
  if (!is.null(obj$stat) && length(obj$stat) != 0) {
    for (i in 1:length(obj$stat)) {
      if (length(obj$stat[[i]]) == length(left.tag))
        obj$stat[[i]] <- obj$stat[[i]][left.tag]
    }
  }
  return (obj)
}

# calcAUCValue
# calculate AUC value with TCC class tcc.
calcAUCValue <- function(tcc) {
  if (is.null(tcc$simulation$trueDE) || length(tcc$simulation$trueDE) == 0)
    stop("\nTCC::ERROR: No true positive annotations about differential expression genes.\n ")
  if (is.null(tcc$stat$rank) || length(tcc$stat$rank) == 0)
    stop("\nTCC::ERROR: There are no rank informations in TCC tcc. It need run TCC.estimateDE().\n")
  return(AUC(rocdemo.sca(truth = as.numeric(tcc$simulation$trueDE != 0), data = - tcc$stat$rank)))
}

# getNormalizedData
# normalize count data with the normalization factors in TCC and return it.
getNormalizedData <- function(tcc) {
  return (tcc$getNormalizedCount())
}


setMethod(
  f = "names",
  signature(x = "TCC"),
  definition = function(x) {
    return (c("count", "names", "group", "norm.factors", 
              "DEGES", "stat", "estimatedDEG", "simulation"))
  }
)
setMethod(
  f = "length",
  signature(x = "TCC"),
  definition = function(x) {
    return (8)
  }
)
setMethod(
  f = "[",
  signature(x = "TCC"),
  definition = function(x, i) {
    if (is.numeric(i)) {
      if (i == 1)
        return (x$count)
      else if (i == 2)
         return (x$names)
      else if (i == 3)
        return (x$group)
      else if (i == 4)
        return (x$norm.factors)
      else if (i == 5)
        return (x$DEES)
      else if (i == 6)
        return (x$stat)
      else if (i == 7)
        return (x$estimatedDEG)
      else if (i == 8)
        return (x$simulation)
      else
        return (NULL)
    } else  {
      if (i == "count")
        return (x$count)
      else if (i == "names")
         return (x$names)
      else if (i == "group")
        return (x$group)
      else if (i == "norm.factors")
        return (x$norm.factors)
      else if (i == "DEGES")
        return (x$DEGES)
      else if (i == "stat")
        return (x$stat)
      else if (i == "estimatedDEG")
        return (x$estimatedDEG)
      else if (i == "simulation")
        return (x$simulation)
      else
        return (NULL)
    }
  }
)

setMethod(
  f = "show",
  signature(object = "TCC"),
  definition = function(object) {
    # Counts.
    cat("Count:\n")
    print(head(object$count))
    cat("\n")
    # Conditions and Annotations.
    df <- data.frame(
      norm.factors = object$norm.factors,
      lib.sizes = object$norm.factors * colSums(object$count)
    )
    rownames(df) <- colnames(object$count)
    df <- cbind(object$group, df)
    cat("Sample:\n")
    print(df)
    cat("\n")
    # Normalized results.
    if (object$private$normalized) {
      cat("DEGES:\n")
      cat(paste("   Pipeline       : ", object$DEGES$protocol, "\n", sep = ""))
      cat(paste("   Execution time : ", sprintf("%.1f", object$DEGES$execution.time[3]),
                " sec\n", sep = ""))
      cat(paste("   Threshold type : ", object$DEGES$threshold$type, 
                " < ", sprintf("%.2f", object$DEGES$threshold$input),"\n",
                "   Potential PDEG : ", sprintf("%.2f", object$DEGES$threshold$PDEG), "\n\n", sep = ""))
    }
    # Esimated results.
    if (object$private$estimated) {
      df <- getResult(object)
      cat("Results:\n")
      print(head(df))
      cat("\n")
    }
  }
)
