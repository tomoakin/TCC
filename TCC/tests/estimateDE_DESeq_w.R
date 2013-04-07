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
cds <- newCountDataSet(tcc$count, group)
sizeFactors(cds) <- tcc$norm.factors * colSums(tcc$count)
cds <- estimateDispersions(cds)
reduced.model <- fitNbinomGLMs(cds, fit0)
full.model <- fitNbinomGLMs(cds, fit1)
p.value <- nbinomGLMTest(full.model, reduced.model)
p.value[is.na(p.value)] <- 1
q.value <- p.adjust(p.value, method = "BH")
tmp <- cbind(p.value, q.value)
rownames(tmp) <- tcc$gene_id
result <- tmp[order(p.value), ]
head(result)
sum(q.value < 0.1)
sum(q.value < 0.2)
