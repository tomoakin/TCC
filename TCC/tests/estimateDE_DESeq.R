library(TCC)
data(hypoData_mg)
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
head(result)
table(tcc$estimatedDEG)
