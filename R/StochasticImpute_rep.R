#create a stochastic impuation function to easily create results over different missing rates
StochasticImpute_rep <- function(data, missing_rate)
{
  X.raw=data #create dataframe that is just the marker data
  X=X.raw#create a copy of the dataframe to work on
  #Set variables to be used to randomize missing values
  n=nrow(X.raw) #rows of data
  m=ncol(X.raw) #columns of data
  dp=m*n #dp gives the total # of data points - use for % calculation later
  uv=runif(dp) #uniform variable
  #We are now ready to set each missing rate
  mr1=missing_rate# missing rate (how much percentage of data is NA)
  #Converts each uniform vector into a series of TRUE/FALSE where each TRUE is a missing variable
  missing1=uv<mr1
  #Returns the true rate of missing generated to ensure function worked
  length(grep("TRUE",missing1))/length(missing1)
  #Creates a matrix of true/false values corresponding to dataframe X in layout
  index.m=matrix(missing1,n,m)#Format indicator as matrix
  #Changes "missing" (TRUE) values to read "NA" and puts them in matrix/dataframe X
  X[index.m]=NA
  n=nrow(X)
  m=ncol(X)
  fn=colSums(X,na.rm=T) # sum of genotypes for all individuals
  fc=colSums(floor(X/3+1),na.rm=T) #count number of non missing individuals
  fa=fn/(2*fc) #Frequency of allele "2"
  for(i in 1:m){
    index.a=runif(n)<fa[i]#if above requency and na include in index 2
    index.na=is.na(X[,i])
    index.m2=index.a&index.na
    index.m0=!index.a&index.na #if not in index 2 include in index 0
    X[index.m2,i]=2 #impute 2
    X[index.m0,i]=0 #impute 1
  }
  accuracy.r=cor(X.raw[index.m],X[index.m], use = 'complete.obs')#correlation accuracy
  index.match=X.raw==X #match index
  index.mm=index.match&index.m
  accuracy.m=length(X[index.mm])/length(X[index.m]) #proportion of match accuracy
  accsi=c(accuracy.r, accuracy.m)
  return(accsi)
}
