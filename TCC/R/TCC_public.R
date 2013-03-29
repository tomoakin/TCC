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
                        FDR=NULL, floorPDEG=NULL, samplesize=10000, processors=NULL) {
      ex.time <- proc.time()
      obj <- tcc$copy()
      obj$calcNormFactors(norm.method=norm.method, test.method=test.method, iteration=iteration,
                      FDR=FDR, floorPDEG=floorPDEG, samplesize=samplesize, processors=processors)
      obj$stat$execution.time <- proc.time() - ex.time
      return(obj)
    }
)

# estimateDE
# the method is for estimating DEGs.
estimateDE <- function(tcc, test.method=NULL, FDR=NULL, samplesize=10000, processors=NULL) {
  obj <- tcc$copy()
  obj$estimateDE(test.method=test.method, FDR=FDR, samplesize=samplesize, processors=processors)
  return(obj)
}

# plot
# plot MA-plot with TCC class tcc.
plot.TCC <- function(x, FDR=NULL, median.lines = FALSE, floor=0, main=NULL, 
                    xlab = expression(A == (log[2] * G2 + log[2] * G1 ) / 2),
                    ylab = expression(M == log[2] * G2 - log[2] * G1),
                    xlim = NULL, ylim = NULL, cex = 0.3, pch = 19, col = NULL, ...) {
      invisible(x$plotMA(FDR=FDR, median.lines=median.lines, floor=floor, main=main, xlab=xlab, ylab=ylab,
               xlim=xlim, ylim=ylim, cex=cex, pch=pch, col=col, ...))
}

# getResult
# get p-value, FDR or the axes of MA-plot as data.frame.
getResult <- function(tcc, sort = FALSE, floor = 0) {
  if (length(table(tcc$group)) != 2)
    stop("\nTCC::EEROR: This version doesn't support when the group more than two.\n")
  if (length(tcc$stat) == 0)
    stop("\nTCC::ERROR: There are no statistics in stat fields of TCC class tcc. Execute TCC.estiamteDE for calculating them.\n")
  count.normed <- tcc$getNormalizedCount()
  mean.exp <- matrix(0, ncol=length(tcc$group), nrow=nrow(tcc$count))
  for (g in 1:length(tcc$group)) {
    mean.exp[, g] <- rowMeans(as.matrix(count.normed[, tcc$group == g]))
  }
  ma.axes <- tcc$.getMACoordinates(mean.exp[, 1], mean.exp[, 2], floor)
  result.df <- data.frame(
    id = rownames(tcc$count),
    a.value = ma.axes$a.value,
    m.value = ma.axes$m.value,
    p.value = tcc$stat$p.value,
    q.value = tcc$stat$q.value,
    rank = tcc$stat$rank,
    estimatedDEG = tcc$estimatedDEG
  )
  if (sort)
    result.df <- result.df[order(result.df$rank), ]
  return (result.df)
}

# filterData
# remove the low count data.
filterLowCountGenes <- function(tcc, low.count = 0) {
  obj <- tcc$copy()
  filters <- matrix(0, ncol=length(obj$group), nrow=nrow(obj$count)) 
  replicates=table(obj$group)
  for (i in 1:length(replicates)) {
    if (replicates[i] == 1) {
      filters[, i] <- as.numeric(obj$count[, (obj$group == i)] <= low.count)
    } else {
      filters[, i] <- as.numeric(rowSums(obj$count[, (obj$group == i)]) <= low.count)
    }
  }
  left.tag <- as.logical(rowSums(filters) != length(replicates))
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
