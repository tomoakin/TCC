library(TCC)
data(hypoData)
df<-data.frame(row.names = paste('a', rownames(hypoData), sep=""),
  A1 = hypoData[,1], A2 = hypoData[,2], A3 = hypoData[,3],
  B1 = hypoData[,4], B2 = hypoData[,5], B3 = hypoData[,6])
m<-cbind(df[[1]], df[[2]], df[[3]],df[[4]],df[[5]],df[[6]])
tccdata=new("TCC", m, c(1, 1, 1, 2, 2, 2))
head(tccdata$gene_id)
