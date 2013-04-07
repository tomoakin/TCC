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
tcc <- estimateDE(tcc, test.method = "edger",
                  FDR = 0.1, design = design, coef = coef)
result <- getResult(tcc, sort = TRUE)
head(result)
table(tcc$estimatedDEG)

