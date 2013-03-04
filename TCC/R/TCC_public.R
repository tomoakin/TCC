# getSimulationData
# sample the simulation data under NB model.
generateSimulationData <- function(Ngene=10000, PDEG=0.20, DEG.assign=c(0.9, 0.1),
                                   DEG.model="uniform", DEG.foldchange=NULL,
                                   group=c(3, 3)) {
# The method is for generating simulation data.
# 1) Make super dispersion from arab data for generating simulation data.
# 2) Adjust disersion$mean for resampling.
#    If "uniform" model of DEG.model, then time foldchange to dispersion$mean.
#    If "gamma" model of DEG.model, then time one to dispersion$mean.
# 3) Generate simulation data under NB dispersion with dispersion$mean.
# 4) Adjust simulation data.
#    If "uniform" model of DEG.model, then times one to all simulation data.
#    If "gamma" model of DEG.model, then times foldchange calculated from DEG.gamma parameters.
# 5) Return the simulation data as matrix object.

  # Prepare and adjust default paramaters.
  max.len <- max(length(DEG.assign), length(group), length(DEG.foldchange))
  if (length(group) != max.len) {
    g <- rep(group, length = max.len)
  } else {
    g <- group
  }
  if (length(DEG.assign) != max.len) {
    def.num <- max.len - length(DEG.assign)
    DEG.assign <- c(DEG.assign[1:(length(DEG.assign) - 1)], 
      rep(DEG.assign[length(DEG.assign)] / (def.num + 1), times=def.num + 1))
  }
  if (is.null(DEG.foldchange)) {
    if (DEG.model == "uniform")
      DEG.foldchange <- list(c(4))
    if (DEG.model == "gamma")
      DEG.foldchange <- list(c(1.2, 2.0, 0.5))
  }
  if (DEG.model == "uniform") {
    for (i in 1:length(DEG.foldchange)) {
      if (length(DEG.foldchange[[i]]) != 1)
        message ("TCC::INFO: DEG.foldchange has three element in the vectors, only the first element is used for fixed foldchange.")
    }
  } else if (DEG.model == "gamma") {
    for (i in 1:length(DEG.foldchange)) {
      if (length(DEG.foldchange[[i]]) != 3)
        stop ("\nTCC::ERROR: It need three elements in each vectors when the DEG.mode is specified to gamma.\n")
    }
  }
  if (length(DEG.foldchange) != max.len) {
    DEG.foldchange <- rep(DEG.foldchange, length=max.len)
  }
    DEG.foldchange <- rep(DEG.foldchange, length = max.len)
  if (sum(DEG.assign) > 1)
    stop("TCC::ERROR: The total value of DEG.assign must less than one.\n") 
  message("TCC::INFO: Generating simulation data under NB distribution ...")
  message(paste("TCC::INFO: (genesizes   : ", paste(Ngene, collapse=", "), ")"))
  message(paste("TCC::INFO: (groups      : ", paste(g, collapse=", "), ")"))
  message(paste("TCC::INFO: (foldhcange distribution : ", DEG.model, ")"))
  message(paste("TCC::INFO: (PDEG        : ", paste(PDEG * DEG.assign, collapse=", "), ")"))

  # 1) Prepare the super population for sampling.
  arab <- NULL
  rm(arab)
  data(arab)
  rpm.a <- sweep(arab[, 1:3], 2, median(colSums(arab[, 1:3])) / colSums(arab[, 1:3]), "*")
  rpm.b <- sweep(arab[, 4:6], 2, median(colSums(arab[, 4:6])) / colSums(arab[, 4:6]), "*")
  rpm.a <- rpm.a[apply(rpm.a, 1, var) > 0, ]
  rpm.b <- rpm.b[apply(rpm.b, 1, var) > 0, ]
  mean.ab <- c(apply(rpm.a, 1, mean), apply(rpm.b, 1, mean))
  var.ab  <- c(apply(rpm.a, 1, var), apply(rpm.b, 1, var))
  dispersion <- (var.ab - mean.ab) / (mean.ab * mean.ab)
  population <- data.frame(mean = mean.ab, disp = dispersion)
  population <- population[population$disp > 0, ]
  resampling.vector <- sample(1:nrow(population), Ngene, replace = TRUE)
  population <- population[resampling.vector, ]  # super dispersion

  # 2) Make foldchagen-matrix for sampling count data.
  fc.matrix <- matrix(1, ncol=sum(g), nrow=Ngene)
  DEG.index <- rep(0, length = nrow(population))              # The DEGs position.
  reps <- rep(1:length(g), times=g)
  if (DEG.model == "uniform") {
    DEG.index[1:round(Ngene * PDEG)] <- 
      rep(1:length(DEG.assign), times = round(Ngene * PDEG * DEG.assign))
    for (i in 1:length(reps)) {
      fc.matrix[, i] <- rep(1, length=Ngene)
      fc.matrix[(DEG.index == reps[i]), i] <- DEG.foldchange[[reps[i]]][1]
    }
  }

  # 3) Sample simulation data from NB dispersion.
  count <- matrix(0, ncol = sum(g), nrow = nrow(population))
  for (i in 1:length(reps)) {
    count[, i] <- rnbinom(n = Ngene, 
      mu = fc.matrix[, i] * population$mean, 
      size = 1 / population$disp)
  }

  # 4) Adjust count data with DEG.gamma paramaters only for "gamma" model.
  if (DEG.model == "gamma") {
    count.means <- matrix(0, ncol=length(g), nrow=Ngene)
    for (i in 1:length(g)) {
      if (is.null(ncol(count[, (reps == i)]))) {
        count.means[, i] <- count[, (reps == i)]
      } else {
        count.means[, i] <- rowMeans(count[, (reps == i)])
      }
    }
    col.idx <- 1
    for (i in 1:length(g)) {
      deg.num <- round(Ngene * PDEG * DEG.assign[i])
      if (is.null(ncol(count.means[, -i]))) {
        deg.candidate <- (count.means[, i] > count.means[, -i])
      } else {
        deg.candidate <- (count.means[, i] > apply(count.means[, -i], 1, max))
      }
      DEG.index[(deg.candidate & cumsum(deg.candidate) <= deg.num)] <- i
      for (j in 1:g[i]) {
        fc.matrix[(DEG.index == i), col.idx] <- 
          DEG.foldchange[[i]][1] + rgamma(sum(DEG.index == i), shape=DEG.foldchange[[i]][2], scale=DEG.foldchange[[i]][3])
        count[(DEG.index == i), col.idx] <- 
          count[(DEG.index == i), col.idx] * fc.matrix[(DEG.index == i), col.idx]
        col.idx <- col.idx + 1
      }
    }
    # sort by DEG.index .
    DEG.index[(DEG.index == 0)] <- 100
    count <- count[order(DEG.index), ]
    fc.matrix <- fc.matrix[order(DEG.index), ]
    DEG.index <- DEG.index[order(DEG.index)]
    DEG.index[(DEG.index == 100)] <- 0
  }

  # save the annotations for generating simulation data to TCC object.
  colnames(count) <- paste("G", rep(1:length(g), times=g), "_rep", sequence(g), sep="")
  rownames(count) <- paste("gene", 1:nrow(count), sep="_") 
  tcc <- new("TCC", count, group)
  tcc$simulation$trueDEG <- DEG.index
  tcc$simulation$DEG.foldchange <- fc.matrix
  tcc$simulation$PDEG <- PDEG * DEG.assign
  tcc$simulation$group <- g
  tcc$replicates <- rep(1:length(g), times=g)
  tcc$private$simulation <- TRUE
  tcc$private$estimated <- FALSE
  return(tcc)
}


# plotSimulationMap
# plot heat map with simulation conditions.
plotFCPseudocolor <- function(tcc, main="",
  xlab="samples", ylab="genes") {
  if (is.null(tcc$simulation$trueDEG) || length(tcc$simulation$trueDEG) == 0)
    message("\nTCC::ERROR: There is no annotations about simulation data.\n")
  # make matrix data for plot heatmap of foldchange.
  d <- tcc$simulation$DEG.foldchange
  # prepare layout.
  layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
  #colorRamp <- rgb(seq(0,1,length=256), seq(0,1,length=256), seq(1,0,length=256))
  maxlevel <- round(max(tcc$simulation$DEG.foldchange))
  colorRamp <- rgb(seq(1, 1, length=maxlevel), seq(1, 0, length=maxlevel), seq(1, 1, length=maxlevel))
  colorLevels <- seq(1, maxlevel, length=length(colorRamp))
  par(mar=c(5.5,4.5,2.5,2))
  image(1:ncol(d), 1:nrow(d), t(d[rev(1:nrow(d)), ]), col=colorRamp,
    ylab=ylab, xlab="", main=main, axes=FALSE, zlim=c(1, max(tcc$simulation$DEG.foldchange)))
  title(xlab=xlab, line=4)
  axis(1, at=1:ncol(d), labels=paste("rep", sequence(tcc$simulation$group), sep=""), cex.axis=0.7, line=0)
  axis(1, at=cumsum(tcc$simulation$group) - tcc$simulation$group + 1,
    labels=paste("Group", c(1:length(tcc$simulation$group)), sep=" "), cex.axis=0.7, line=1, lty=0)
  y.axis <- c(1, cumsum(nrow(tcc$count) * tcc$simulation$PDEG), nrow(tcc$count) - 0.5)
  y.labels <- c(1, cumsum(nrow(tcc$count) * tcc$simulation$PDEG), nrow(tcc$count))
  axis(2, at=nrow(tcc$count) - y.axis, labels=y.labels, cex.axis=0.7, las=1)
  box()
  # colorbar.
  par(mar = c(5.5, 2.5, 2.5, 2))
  image(1, colorLevels, matrix(colorLevels, ncol=length(colorLevels), nrow=1),
    col = colorRamp, xlab="", ylab="", xaxt="n")
  box()
  layout(1)
}



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
                    xlim = NULL, ylim = NULL, cex = 0.3, pch = 19, col = NULL, col.tag = NULL,...) {
      invisible(x$plotMA(FDR=FDR, median.lines=median.lines, floor=floor, main=main, xlab=xlab, ylab=ylab,
               xlim=xlim, ylim=ylim, cex=cex, pch=pch, col=col, col.tag=col.tag,...))
}

# getResult
# get p-value, FDR or the axes of MA-plot as data.frame.
getResult <- function(tcc, sort = FALSE, floor = 0) {
  if (length(tcc$group) != 2)
    stop("\nTCC::EEROR: This version doesn't support when the group more than two.\n")
  if (length(tcc$stat) == 0)
    stop("\nTCC::ERROR: There are no statistics in stat fields of TCC class tcc. Execute TCC.estiamteDE for calculating them.\n")
  count.normed <- tcc$getNormalizedCount()
  mean.exp <- matrix(0, ncol=length(tcc$group), nrow=nrow(tcc$count))
  for (g in 1:length(tcc$group)) {
    mean.exp[, g] <- rowMeans(as.matrix(count.normed[, tcc$replicates == g]))
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
  for (i in 1:length(obj$group)) {
    if (obj$group[i] == 1) {
      filters[, i] <- as.numeric(obj$count[, (obj$replicates == i)] <= low.count)
    } else {
      filters[, i] <- as.numeric(rowSums(obj$count[, (obj$replicates == i)]) <= low.count)
    }
  }
  left.tag <- as.logical(rowSums(filters) != length(obj$group))
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
    nm <- c("count", "names", "group", "replicates", "norm.factors")
    if (!is.null(x$stat))
      nm <- c(nm, "stat")
    if (!is.null(x$estimatedDEG))
      nm <- c(nm, "estimateDEG")
    if (!is.null(x$simulation))
      nm <- c(nm, "simulation")
  }
)

setMethod(
  f = "show",
  signature(object = "TCC"),
  definition = function(object) {
    if (object$private$estimated) {
      df <- getResult(object)
      cat("Analyzed results:\n")
      print(head(df))
      cat("\n")
    } else {
      cat("Count:\n")
      print(head(object$count))
      cat("\n")
      df <- data.frame(
        replicates = object$replicates,
        norm.factors = object$norm.factors
      )
      rownames(df) <- colnames(object$count)
      cat("Condition and normalization factors:\n")
      print(df)
      cat("\n")
    }
  }
)
