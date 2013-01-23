library(TCC)
data(hypoData)
df<-data.frame(names = hypoData[,0],
  A1 = hypoData[,1], A2 = hypoData[,2], A3 = hypoData[,3],
  B1 = hypoData[,4], B2 = hypoData[,5], B3 = hypoData[,6])
tccdata=new("TCC",df,c(3,3))
