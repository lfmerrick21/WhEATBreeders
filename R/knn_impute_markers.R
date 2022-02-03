#create a KNN impuation function to easily create results over different missing rates
knn_impute_markers <- function(data, missing_rate)
{
  X.raw<-data#create dataframe that is just the marker data
  X <- X.raw#create a copy of the dataframe to work on
  #Set variables to be used to randomize missing values
  mr=missing_rate # missing rate (how much percentage of data is NA)
  n=nrow(X)
  m=ncol(X)
  dp=m*n ##dp gives the total # of data points - use for % calculation later
  #Use KNN to impute similar to problem 3
  #Converts each uniform vector into a series of TRUE/FALSE where each TRUE is a missing variable
  uv=runif(dp)
  missing=uv<mr
  index.m=matrix(missing,n,m)#Format indicator as matrix
  X[index.m]=NA #insert NA into index
  x=impute::impute.knn(as.matrix(t(X)))#transpose and impute
  accuracy.r=cor(X.raw[index.m],t(x$data)[index.m], use = 'complete.obs')
  index.match=X.raw==t(x$data)
  index.mm=index.match&index.m
  accuracy.m=length(X[index.mm])/length(X[index.m])
  accknn=c(accuracy.r, accuracy.m)
  return(accknn)
}
