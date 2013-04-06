selectMethodByGroup <- function(group){
  if ((ncol(group) == 1) && (min(as.numeric(table(group))) == 1)) 
    "deseq"
  else 
    "edger"
}
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
                                cl = NULL){
  ex.time <- proc.time()
  if (is.null(norm.method)) 
      norm.method <- selectMethodByGroup(group)
  if (is.null(test.method))
      test.method <- selectMethodByGroup(group)
  if (norm.method == "edger")
    norm.method <- "tmm" 
  if (test.method != "bayseq"){
    if (is.null(FDR)) FDR <- 0.1
    if (is.null(floorPDEG)) floorPDEG <- 0.05
  }
  if (iteration) {
    if (is.logical(iteration))
      iteration <- 1
    if (iteration < 0 && 100 < iteration)
      stop("TCC::ERROR: The iteration must be given within the limits of from '0' to '100'.")
    message(paste("TCC::INFO: Calculating normalization factors using DEGES"))
    message(paste("TCC::INFO: (iDEGES pipeline :", norm.method, "- [", test.method, "-", norm.method, "] X", iteration, ")"))
    DEGES$protocol <<- paste(norm.method, "- [", test.method, "-", norm.method, "] X", iteration)
  } else {
    message(paste("TCC::INFO: Calculating normalization factors using", norm.method, "..."))
    DEGES$protocol <<- norm.method
  }
  # DEGES strategy STEP 1. (First normalization)
  norm.factors <<- switch(norm.method,
    "tmm" = .self$.normByTmm(count),
    "deseq" = .self$.normByDeseq(count),
    stop(paste("\nTCC::ERROR: The normalization method of ", norm.method, " doesn't supported.\n"))
  )
  message("TCC::INFO: First normalization done.")
  norm.factors <<- norm.factors / mean(norm.factors)
  DEGES$threshold <<- list(type = "Unused", input = 0, PDEG = 0)
  #  if DEGES not NULL then start DEGES strategy.
  if (iteration > 0) {
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
      # DEGES strategy STEP 3. (Second normalization)
      norm.factors <<- switch(norm.method,
        "tmm" = .self$.normByTmm(count.ndeg),
        "deseq" = .self$.normByDeseq(count.ndeg)
      )
      norm.factors <<- norm.factors * colSums(count.ndeg) / colSums(count)
      norm.factors <<- norm.factors / mean(norm.factors)
    }
    DEGES$potentialDEG <<- deg.flg
  }
  message("TCC::INFO: Done.")
  DEGES$execution.time <<- proc.time() - ex.time
  private$normalized <<- TRUE
})
