library(TCC)
data(hypoData)
set.seed(1834016)

df<-data.frame(row.names = paste('a', rownames(hypoData), sep=""),
  A1 = hypoData[,1], A2 = hypoData[,2], A3 = hypoData[,3],
  B1 = hypoData[,4], B2 = hypoData[,5], B3 = hypoData[,6])
head(df)
tccdata=new("TCC",df,c(1, 1, 1, 2, 2, 2))
head(tccdata$gene_id)

df<-data.frame(row.names = paste('a', rownames(hypoData), sep=""),
  A1 = hypoData[,1], A2 = hypoData[,2], A3 = hypoData[,3],
  B1 = hypoData[,4], B2 = hypoData[,5], B3 = hypoData[,6])
m<-cbind(df[[1]], df[[2]], df[[3]],df[[4]],df[[5]],df[[6]])
head(m)
tccdata=new("TCC", m, c(1, 1, 1, 2, 2, 2))
head(tccdata$gene_id)

group <- c(1, 1, 1, 2, 2, 2)
tcc <- new("TCC", hypoData, group)
cat("tcc$count: ")
cat(dim(tcc$count))
cat("\n")
tccf <- filterLowCountGenes(tcc)
cat("dim(tcc$count): ")
cat(dim(tccf$count))
cat("\n")
cat("dim(hypoData): ")
cat(dim(hypoData))
cat("\n")
cat("dim(hypoData[as.logical(rowSums(hypoData)>0),]): ")
cat(dim(hypoData[as.logical(rowSums(hypoData) > 0),]))
cat("\n")

tcc <- simulateReadCounts(Ngene = 100, PDEG = 0.2, DEG.assign = c(0.9, 0.1))
#  group <- c(3, 3)
#  tcc <- new("TCC", hypoData, c(3,3))
show(tcc)
sub_tcc <- subset(tcc,1:10*10)
show(sub_tcc)
sub2_tcc <- tcc[1:10*10-1]
show(sub2_tcc)
tcc$calcNormFactors(norm.method = "tmm", test.method = "bayseq")
tcc$estimateDE(test.method = "edger", FDR = 0.1)

sub_tcc <- subset(tcc,1:10*10)
show(sub_tcc)

show(tcc[c("gene_13", "gene_17", "gene_23")])
result <- getResult(sub_tcc, sort = TRUE)
print(head(result))

tcc <- simulateReadCounts(Ngene = 100, PDEG = 0.2, DEG.assign = c(0.9, 0.1))
tcc <- calcNormFactors(tcc)
tcc <- estimateDE(tcc, test.method = "bayseq", FDR = 0.1, samplesize = 100)
result <- getResult(tcc, sort = TRUE)
print(head(result))
table(tcc$estimatedDEG) 
png("plot4b.png", 600, 500)
plot(tcc)
dev.off()

tcc <- simulateReadCounts(Ngene = 100, PDEG = 0.2, DEG.assign = c(0.9, 0.1))
tcc <- calcNormFactors(tcc)
tcc <- estimateDE(tcc, test.method = "edger", FDR = 0.1)
result <- getResult(tcc, sort = TRUE)
print(head(result))
table(tcc$estimatedDEG) 
png("plot4.png", 600, 500)
plot(tcc)
dev.off()

data(hypoData_mg)
group <- c(1, 1, 1, 2, 2, 2, 3, 3, 3)
tcc <- new("TCC", hypoData_mg, group)
design <- model.matrix(~ as.factor(group))
coef <- 2:length(unique(group))
tcc <- calcNormFactors(tcc, norm.method = "tmm", test.method = "edger",
                        iteration = 1, design = design, coef = coef)
#tcc$norm.factors
## this will differ depending on edgeR versions

tcc <- simulateReadCounts(Ngene = 100, PDEG = 0.2, DEG.assign = c(0.9, 0.1),replicate=c(1,1))
tcc <- calcNormFactors(tcc, norm.method = "deseq", test.method = "deseq",
                          FDR = 0.1, floorPDEG = 0.05)
cat("tcc$norm.factors: ")
cat(tcc$norm.factors)
cat("\n")
cat("tcc$stat$execution.time: ")
cat(tcc$stat$execution.time)
cat("\n")

tcc2 <- simulateReadCounts(Ngene = 100, PDEG = 0.2, DEG.assign = c(0.9, 0.1))
tcc <- calcNormFactors(tcc2, norm.method = "deseq", test.method = "deseq",
                        FDR = 0.1, floorPDEG = 0.05)
cat("tcc$norm.factors: ")
cat(tcc$norm.factors)
cat("\n")
cat("tcc$stat$execution.time: ")
cat(tcc$stat$execution.time)
cat("\n")

tcc <- calcNormFactors(tcc2, iteration = 3)
cat("tcc$norm.factors: ")
cat(tcc$norm.factors)
cat("\n")
cat("tcc$stat$execution.time: ")
cat(tcc$stat$execution.time)
cat("\n")

tcc <- calcNormFactors(tcc2, norm.method = "tmm", test.method = "edger")
cat("tcc$norm.factors: ")
cat(tcc$norm.factors)
cat("\n")
cat("tcc$stat$execution.time: ")
cat(tcc$stat$execution.time)
cat("\n")

tcc <- calcNormFactors(tcc2, norm.method = "tmm", test.method = "bayseq")
cat("tcc$norm.factors: ")
cat(tcc$norm.factors)
cat("\n")
cat("tcc$stat$execution.time: ")
cat(tcc$stat$execution.time)
cat("\n")

group <- c(1, 1, 1, 2, 2, 2, 3, 3, 3)
tcc <- new("TCC", hypoData_mg, group)
###  Normalization  ###
design <- model.matrix(~ as.factor(group))
coef <- 2:length(unique(group))
tcc <- calcNormFactors(tcc, norm.method = "tmm", test.method = "edger",
                        iteration = 1, design = design, coef = coef)
###  DE analysis  ###
fit1 <- count ~ condition
fit0 <- count ~ 1
tcc <- estimateDE(tcc, test.method = "deseq",
                  FDR = 0.1, fit0 = fit0, fit1 = fit1)
result <- getResult(tcc, sort = TRUE)
#head(result)
## this will differ depending on edgeR versions
#table(tcc$estimatedDEG)

group <- c(1, 1, 1, 2, 2, 2, 3, 3, 3)
tcc <- new("TCC", hypoData_mg, group)
###  Normalization  ###
design <- model.matrix(~ as.factor(group))
coef <- 2:length(unique(group))
tcc <- calcNormFactors(tcc, norm.method = "tmm", test.method = "edger",
                        iteration = 1, design = design, coef = coef)

group <- c(1, 1, 1, 2, 2, 2, 3, 3, 3)
tcc <- new("TCC", hypoData_mg, group)
###  Normalization  ###
design <- model.matrix(~ as.factor(group))
coef <- 2:length(unique(group))
tcc <- calcNormFactors(tcc, norm.method = "tmm", test.method = "edger",
                        iteration = 1, design = design, coef = coef)
###  DE analysis  ###
tcc <- estimateDE(tcc, test.method = "edger",
                  FDR = 0.1, design = design, coef = coef)
result <- getResult(tcc, sort = TRUE)
#head(result)
## this will differ depending on edgeR versions
#table(tcc$estimatedDEG)

tcc <- simulateReadCounts(Ngene = 100, PDEG = 0.3,
                       DEG.assign = c(0.6, 0.2, 0.2),
                       DEG.foldchange = c(3, 10, 6),
                       replicates = c(2, 4, 3))
cat("dim(tcc$count): ")
cat(dim(tcc$count))
cat("\n")
cat("tcc$group$group: ")
cat(tcc$group$group)
cat("\n")
cat("tcc$count:\n")
print(head(tcc$count))
plotFCPseudocolor(tcc)
plot(tcc)

tcc <- simulateReadCounts(Ngene = 200, PDEG = 0.30,
                              DEG.assign = c(0.85, 0.15),
                              DEG.foldchange = list(min = c(1.2, 1.2),
                                                    shape = c(2.0, 2.0),
                                                    scale = c(0.5, 0.5)),
                              replicates = c(2, 2))
tcc <- calcNormFactors(tcc, norm.method = "tmm", test.method = "edger",
                       iteration = 3, FDR = 0.1, floorPDEG = 0.05)
plot(tcc, median.lines = TRUE)

