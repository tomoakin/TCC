# getSimulationData
# sample the simulation data under NB model.
generateSimulationData <- function(Ngene=10000, PDEG=0.20, DEG.assign=c(0.9, 0.1),
                                   DEG.model="uniform", DEG.foldchange=NULL,
                                   replicates=c(3, 3)) {
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
  if (class(DEG.foldchange) == "list") {
    max.len <- max(length(DEG.assign), length(replicates), length(DEG.foldchange[[1]]))
  } else {
    max.len <- max(length(DEG.assign), length(replicates), length(DEG.foldchange))
  }
  if (length(replicates) != max.len) {
    g <- rep(replicates, length = max.len)
  } else {
    g <- replicates
  }
  if (length(DEG.assign) != max.len) {
    def.num <- max.len - length(DEG.assign)
    DEG.assign <- c(DEG.assign[1:(length(DEG.assign) - 1)], 
      rep(DEG.assign[length(DEG.assign)] / (def.num + 1), times=def.num + 1))
  }
  if (is.null(DEG.foldchange)) {
    if (DEG.model == "uniform")
      DEG.foldchange <- rep(4, length = max.len)
    if (DEG.model == "gamma")
      DEG.foldchange <- lapply(list(1.2, 2.0, 0.5), function(l){rep(l, length = max.len)})
  }
  if (DEG.model == "gamma" && length(DEG.foldchange) != 3)
    stop ("\nTCC::ERROR: It need a list object contained three vectors when the DEG.mode is specified to gamma.\n")

  if (sum(DEG.assign) > 1)
    stop("TCC::ERROR: The total value of DEG.assign must less than one.\n") 
  message("TCC::INFO: Generating simulation data under NB distribution ...")
  message(paste("TCC::INFO: (genesizes   : ", paste(Ngene, collapse=", "), ")"))
  message(paste("TCC::INFO: (replicates  : ", paste(g, collapse=", "), ")"))
  message(paste("TCC::INFO: (foldhcange distribution : ", DEG.model, ")"))
  message(paste("TCC::INFO: (PDEG        : ", paste(PDEG * DEG.assign, collapse=", "), ")"))

  # 1) Prepare the super population for sampling.
  group <- rep(1:length(g), times = g)
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
      fc.matrix[(DEG.index == reps[i]), i] <- DEG.foldchange[reps[i]]
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
          DEG.foldchange[[1]][i] + rgamma(sum(DEG.index == i), shape=DEG.foldchange[[2]][i], scale=DEG.foldchange[[3]][i])
        count[(DEG.index == i), col.idx] <- 
          count[(DEG.index == i), col.idx] * fc.matrix[(DEG.index == i), col.idx]
        col.idx <- col.idx + 1
      }
      count <- round(count)
    }
    # sort by DEG.index .
    DEG.index[(DEG.index == 0)] <- 100
    count <- count[order(DEG.index), ]
    fc.matrix <- fc.matrix[order(DEG.index), ]
    DEG.index <- DEG.index[order(DEG.index)]
    DEG.index[(DEG.index == 100)] <- 0
  }
  colnames(count) <- paste("G", rep(1:length(g), times=g), "_rep", sequence(g), sep="")
  rownames(count) <- paste("gene", 1:nrow(count), sep="_") 
  # Adjust column index.
  #count.adjust <- matrix(0, ncol = ncol(count), nrow = nrow(count))
  #fc.matrix.adjust <- matrix(0, ncol = ncol(count), nrow = nrow(count))
  #labels.old <- rep(unique(group), times = replicates)
  #labels <- table(group)
  #for (i in 1:length(labels)) {
  #  count.adjust[, (group == names(labels)[[i]])] <- count[, (labels.old == names(labels)[[i]])]
  #  fc.matrix.adjust[, (group == names(labels)[[i]])] <- fc.matrix[, (labels.old == names(labels)[[i]])]
  #}
  tcc <- new("TCC", count, group)
  tcc$simulation$trueDEG <- DEG.index
  tcc$simulation$DEG.foldchange <- fc.matrix
  tcc$simulation$PDEG <- PDEG * DEG.assign
  tcc$private$simulation.rep <- g
  tcc$private$simulation <- TRUE
  tcc$private$estimated <- FALSE
  return(tcc)
}


# plotSimulationMap
# plot heat map with simulation conditions.
plotFCPseudocolor <- function(tcc, main="", xlab="samples", ylab="genes") {
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
  axis(1, at=1:ncol(d), labels=paste("rep", sequence(tcc$private$simulation.rep), sep=""), cex.axis=0.7, line=0)
  axis(1, at=cumsum(tcc$private$simulation.rep) - tcc$private$simulation.rep + 1,
    labels=paste("Group", c(1:length(tcc$private$simulation.rep)), sep=" "), cex.axis=0.7, line=1, lty=0)
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


