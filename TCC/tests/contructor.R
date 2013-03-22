library(TCC)
data(hypoData)

# The data frame with row names.
df<-data.frame(row.names = paste('a', rownames(hypoData), sep=""),
               A1 = hypoData[, 1], A2 = hypoData[, 2], A3 = hypoData[, 3],
               B1 = hypoData[, 4], B2 = hypoData[, 5], B3 = hypoData[, 6])
tccdata <- new("TCC", df, c(1, 1, 1, 2, 2, 2))
head(tccdata$count)

# The data frame without row names.
df<-data.frame(A1 = hypoData[, 1], A2 = hypoData[, 2], A3 = hypoData[, 3],
               B1 = hypoData[, 4], B2 = hypoData[, 5], B3 = hypoData[, 6])
rownames(df) <- NULL
tccdata <- new("TCC", df, c(1, 1, 1, 2, 2, 2))
head(tccdata$count)

# The matrix with row names.
mt <- hypoData
tccdata <- new("TCC", mt, c(1, 1, 1, 2, 2, 2))
head(tccdata$count)

# The matrix without row names.
mt <- hypoData
rownames(mt) <- NULL
tccdata <- new("TCC", mt, c(1, 1, 1, 2, 2, 2))
head(tccdata$count)

# Multiple groups and mutiple factors
mt <- hypoData
rownames(mt) <- NULL
group <- data.frame(replicates = c(1, 1, 1, 2, 2, 2), 
                    times = c(1, 2, 3, 1, 2, 3))
tccdata <- new("TCC", mt, group)

