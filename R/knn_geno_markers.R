#create a KNN marker impuation function to easily create results over different missing rates
knn_geno_markers <- function(data, missing_rate)
{
  X.rawrev=t(data)#create dataframe that is just the marker data tranposed
  X.rev=X.rawrev#create a copy of the dataframe to work on
  #Sets the constants we will need (which are different than 3!)
  n=nrow(X.rawrev)
  m=ncol(X.rawrev)
  dp=m*n #dp gives the total # of data points - use for % calculation later
  #Use KNN to impute similar to problem 3
  X.rev=X.rawrev
  mr=missing_rate# missing rate (how much percentage of data is NA)
  uv=runif(dp) #Converts each uniform vector into a series of TRUE/FALSE where each TRUE is a missing variable
  missing=uv<mr
  index.m=matrix(missing,n,m)#Format indicator as matrix
  X.rev[index.m]=NA
  x=impute::impute.knn(as.matrix(t(X.rev))) #impute transpose of transposed
  accuracy.r=cor(X.rawrev[index.m],t(x$data)[index.m], use = 'complete.obs')
  index.match=X.rawrev==t(x$data)
  index.mm=index.match&index.m
  accuracy.m=length(X.rev[index.mm])/length(X.rev[index.m])
  accknngeno=c(accuracy.r, accuracy.m)
  return(accknngeno)
}
